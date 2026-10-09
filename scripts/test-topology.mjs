import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { readFileSync, mkdtempSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { execFileSync } from 'node:child_process';
import { createServer } from 'vite';

// Compare full-precision internal states with the unchanged heterogeneous model.
const parityDirectory = mkdtempSync(join(tmpdir(), 'potts-dynamics-'));
try {
  const outputs = ['hybrid', 'topology'].map(model => {
    const root = `backend/model_potts_${model}`;
    const program = join(parityDirectory, model);
    execFileSync('g++', ['-O3', '-std=gnu++17', `-I${root}/include`,
      '-Ibackend/model_potts_hybrid/include', 'tests/potts-dynamics-parity.cpp',
      `${root}/pnet.cpp`, 'backend/model_potts_hybrid/functions.cpp',
      'backend/model_potts_hybrid/rand_gen.cpp', '-o', program]);
    return execFileSync(program, [], { encoding: 'utf8', stdio: ['ignore', 'pipe', 'ignore'] });
  });
  assert.equal(outputs[1], outputs[0], 'Adaptation and both inhibition states match heterogeneous dynamics');
  assert.equal(outputs[0].trim().split('\n').length, 6*301);
} finally { rmSync(parityDirectory, { recursive: true, force: true }); }

const directory = 'src/lib/topologies';
const manifest = JSON.parse(readFileSync(`${directory}/manifest.json`, 'utf8'));
const source = JSON.parse(readFileSync('tests/topology-source.json', 'utf8'));
assert.equal(manifest.topologies.length, 12);
for (const [file, sha] of Object.entries(source.blobs)) {
  const bytes = readFileSync(`${directory}/${file}`);
  assert.equal(createHash('sha1').update(`blob ${bytes.length}\0`).update(bytes).digest('hex'), sha,
    `Frozen upstream blob: ${file}`);
}
for (const item of manifest.topologies) {
  const matrix = readFileSync(`${directory}/${item.file}`, 'utf8').trim().split('\n').map(row => row.split(',').map(Number));
  assert.equal(matrix.length, 49);
  assert(matrix.every(row => row.length === 49 && row.every(value => value === 0 || value === 1)));
  for (let f = 0; f < 49; ++f) {
    assert.equal(matrix[f].reduce((x, y) => x+y), item.degree);
    assert.equal(matrix.reduce((sum, row) => sum+row[f], 0), item.degree);
  }
  let n4 = 0;
  for (let f = 0; f < 49; ++f) for (let g = f+1; g < 49; ++g) {
    const shared = matrix[f].reduce((sum, edge, p) => sum+edge*matrix[g][p], 0);
    n4 += shared*(shared-1)/2;
  }
  assert.equal(n4, item.metrics.n4);
  assert(Math.abs(Math.log1p(n4)-item.metrics.redundancy_log1p) < 1e-12);
  const args = [`${directory}/${item.file}`, 1.1, 1.1, 7, 7, 100, 400, .5, .5,
    .1, 11, .25, 20, .15, 1, 17, 0, 0, 0, 0, 60, 20].map(String);
  const rows = execFileSync('.simulators/topology', args, { encoding: 'utf8', stdio: ['ignore', 'pipe', 'ignore'] })
    .trim().split('\n').map(line => line.trim().split(/\s+/).map(Number));
  assert.equal(rows.length, 6);
  assert(rows.every(row => row.length === 50 && row.every(Number.isFinite)));
  assert.deepEqual(rows.map(row => row[0]), [0, 0, 10, 10, 20, 20]);
}
const fixtures = JSON.parse(readFileSync('tests/topology-reference.json', 'utf8')).cases;
for (const fixture of fixtures) {
  const item = manifest.topologies.find(item => item.id === fixture.topology);
  const output = execFileSync('.simulators/topology', [`${directory}/${item.file}`, ...fixture.args],
    { encoding: 'utf8', stdio: ['ignore', 'pipe', 'ignore'] });
  assert.deepEqual(output.trim().split('\n').map(line => line.trim().split(/\s+/).map(Number)), fixture.rows,
    `Dense corrected-reference numerical regression: ${fixture.name ?? fixture.topology}`);
}
// Isolated regions make the regional cue routing directly observable.
const independent = fixtures.find(item => item.name === 'independent-stored-cues');
const winner = row => row.indexOf(Math.max(...row.slice(1))) - 1;
assert.equal(winner(independent.rows[2]), 3, 'Posterior receives its own stored cue');
assert.equal(winner(independent.rows[3]), 17, 'Frontal receives its own stored cue');
for (const name of ['posterior-random-cue', 'frontal-random-cue', 'both-random-cues']) {
  const fixture = fixtures.find(item => item.name === name);
  for (const [offset, random] of [[2, fixture.args[17] === '1'], [3, fixture.args[18] === '1']]) {
    const maximum = Math.max(...fixture.rows[offset].slice(1));
    assert(random ? maximum < .85 : maximum > .9, `${name}: external cue is not a stored-memory cue`);
  }
}
const randomFixture = fixtures.find(item => item.name === 'both-random-cues');
const ignoredCueArgs = [...randomFixture.args];
ignoredCueArgs[15] = '48'; ignoredCueArgs[16] = '19';
const ignoredCueRows = execFileSync('.simulators/topology', [`${directory}/T00_one_to_one.csv`, ...ignoredCueArgs],
  { encoding: 'utf8', stdio: ['ignore', 'pipe', 'ignore'] }).trim().split('\n').map(row => row.trim().split(/\s+/).map(Number));
assert.deepEqual(ignoredCueRows, randomFixture.rows, 'Random cue ignores the disabled stored-memory sliders');
console.log('PASS: heterogeneous dynamics parity, exact upstream topology blobs, 49×49 invariants, all 12 native runs, eight corrected dense regressions, independent stored and unstored random cues.');

const server = await createServer({ server: { host: '127.0.0.1', port: 0 } });
await server.listen();
const base = `http://127.0.0.1:${server.httpServer.address().port}`;
const parameters = {
  model: 'topology', topology: 'T01', steps: 30,
  posterior: { S: 3, w: 1.1, tau2: 200, lambda: .5, cue: 0, randomCue: false },
  frontal: { S: 3, w: 1.1, tau2: 200, lambda: .5, cue: 0, randomCue: false },
  global: { U: .1, beta: 11, a: .25, tau1: 20, density: 50/256 }
};
const request = body => fetch(`${base}/api/simulate`, {
  method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body)
});
async function simulate(body) {
  const response = await request(body);
  assert.equal(response.status, 200);
  const events = (await response.text()).trim().split('\n').map(JSON.parse);
  assert.equal(events[0].type, 'started');
  assert.equal(events[0].memories, 49);
  assert.equal(events[0].topology, body.topology);
  assert.equal(events[0].posteriorCue, body.posterior.randomCue ? 'random' : body.posterior.cue);
  assert.equal(events[0].frontalCue, body.frontal.randomCue ? 'random' : body.frontal.cue);
  assert.equal(events.at(-1).type, 'done');
  const samples = events.filter(event => event.type === 'sample');
  assert.equal(samples.length, 4);
  assert.equal(samples.at(-1).time, 30);
  assert(samples.every(sample => sample.posterior.length === 49 && sample.frontal.length === 49 &&
    [...sample.posterior, ...sample.frontal].every(Number.isFinite)));
  return { info: events[0], samples };
}
try {
  const [first, second] = await Promise.all([simulate(parameters), simulate(parameters)]);
  assert.notEqual(first.info.patternSeed, second.info.patternSeed);
  assert.notEqual(first.info.runtimeSeed, second.info.runtimeSeed);
  assert.notDeepEqual(first.samples, second.samples, 'Repeated starts refresh both regions');
  await simulate({ ...parameters, topology: 'T00' });
  await simulate({ ...parameters, topology: 'T11', global: { U: .6, beta: 21, a: .4, tau1: 5, density: .25 },
    posterior: { ...parameters.posterior, cue: 48 }, frontal: { ...parameters.frontal, cue: 17 } });
  for (const [randomP, randomF] of [[true, false], [false, true], [true, true]]) {
    await simulate({ ...parameters,
      posterior: { ...parameters.posterior, randomCue: randomP },
      frontal: { ...parameters.frontal, randomCue: randomF } });
  }
  for (const body of [
    { ...parameters, topology: '../../etc/passwd' }, { ...parameters, topology: 'T12' },
    { ...parameters, global: undefined },
    { ...parameters, posterior: { ...parameters.posterior, lambda: 1.01 } },
    { ...parameters, frontal: { ...parameters.frontal, tau2: 100 } },
    ...[['U', -.1], ['beta', 22], ['beta', '11'], ['a', .5], ['tau1', 0], ['density', .3]].map(([key, value]) =>
      ({ ...parameters, global: { ...parameters.global, [key]: value } })),
    ...['posterior', 'frontal'].flatMap(area =>
      [['cue', 49], ['cue', .5], ['cue', -1], ['cue', '0'], ['cue', undefined], ['randomCue', 'true'], ['randomCue', 1], ['randomCue', undefined]]
        .map(([key, value]) => ({ ...parameters, [area]: { ...parameters[area], [key]: value } })))
  ]) assert.equal((await request(body)).status, 400);
  const html = await (await fetch(`${base}/models/potts_topology`)).text();
  assert.match(html, /Many-to-many Potts Network/);
  assert.match(html, /U = 0\.10/);
  assert.match(html, /256 units and 49/);
  assert.match(html, /cₘ = 50 · C\/N = 0\.19531/);
  assert.equal((html.match(/τ₂ = 200/g) ?? []).length, 2);
  assert.doesNotMatch(html, /Choose a frontal adaptation time|Frontal · slow|Posterior · fast/);
  assert.match(html, /Adaptation tracks state activity σ/);
  assert.match(html, /T₃A = 10, T₃B = 100000 and γA = 0\.5/);
  assert.doesNotMatch(html, /tracks 1\.2r|inhibition is disabled|no active effect/);
  assert.equal((html.match(/<option value="T\d\d"/g) ?? []).length, 12);
  assert.match(html, /score-help-2/);
  assert.equal((html.match(/Cue or random/g) ?? []).length, 2);
  assert.match(html, /aria-label="Frontal random cue"/);
  assert.match(html, /aria-label="Posterior random cue"/);
  assert.equal((html.match(/role="tooltip"/g) ?? []).length, 18);
  assert.doesNotMatch(html, /simcode_pilot|Generator order stored|generator order|Cued memory index/);
  assert.doesNotMatch(html, /id="display-order"/);
  const { topologies } = await server.ssrLoadModule('/src/lib/topologies.ts');
  const { render } = await server.ssrLoadModule('svelte/server');
  const { default: TopologyView } = await server.ssrLoadModule('/src/routes/models/potts_topology/TopologyView.svelte');
  for (const topology of topologies) {
    const csv = readFileSync(`${directory}/${topology.file}`, 'utf8').trim().split('\n').map(row => row.split(',').map(Number));
    assert.deepEqual(topology.matrix, csv, `${topology.id}: loading retains exact CSV order`);
    const view = render(TopologyView, { props: { topology } }).body;
    const cells = [...view.matchAll(/<rect x="([\d.]+)" y="([\d.]+)" width="5.5" height="5.5"/g)]
      .map(([, x, y]) => [Number(y), Number(x)]);
    const expected = csv.flatMap((row, f) => row.flatMap((edge, p) => edge ? [[12+f*6, 12+p*6]] : []));
    assert.deepEqual(cells, expected, `${topology.id}: every SVG cell retains its CSV coordinates`);
  }
  console.log('PASS: fresh concurrent seeds, all regional cue modes, API validation, parameter tooltips, and all 12 rendered matrices match CSVs cell for cell.');
} finally { await server.close(); }
