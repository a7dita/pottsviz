import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { execFileSync } from 'node:child_process';
import { createServer } from 'vite';

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
    .3, 11, .25, 20, .15, 1, 17, 0, 60, 20].map(String);
  const rows = execFileSync('.simulators/topology', args, { encoding: 'utf8', stdio: ['ignore', 'pipe', 'ignore'] })
    .trim().split('\n').map(line => line.trim().split(/\s+/).map(Number));
  assert.equal(rows.length, 6);
  assert(rows.every(row => row.length === 50 && row.every(Number.isFinite)));
  assert.deepEqual(rows.map(row => row[0]), [0, 0, 10, 10, 20, 20]);
}
for (const fixture of JSON.parse(readFileSync('tests/topology-reference.json', 'utf8')).cases) {
  const item = manifest.topologies.find(item => item.id === fixture.topology);
  const output = execFileSync('.simulators/topology', [`${directory}/${item.file}`, ...fixture.args],
    { encoding: 'utf8', stdio: ['ignore', 'pipe', 'ignore'] });
  assert.deepEqual(output.trim().split('\n').map(line => line.trim().split(/\s+/).map(Number)), fixture.rows,
    `Single-threaded pilot numerical regression: ${fixture.topology}`);
}
console.log('PASS: exact upstream topology blobs, 49×49 degrees and redundancy scores, all 12 native runs, four pilot numerical regressions.');

const server = await createServer({ server: { host: '127.0.0.1', port: 0 } });
await server.listen();
const base = `http://127.0.0.1:${server.httpServer.address().port}`;
const parameters = {
  model: 'topology', topology: 'T01', steps: 30,
  posterior: { S: 3, w: 1.1, tau2: 100, lambda: .5 },
  frontal: { S: 3, w: 1.1, tau2: 400, lambda: .5 },
  global: { U: .3, beta: 11, a: .25, tau1: 20, density: .15, cue: 0 }
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
  await simulate({ ...parameters, topology: 'T11', global: { U: .6, beta: 21, a: .4, tau1: 5, density: .25, cue: 48 } });
  for (const body of [
    { ...parameters, topology: '../../etc/passwd' }, { ...parameters, topology: 'T12' },
    { ...parameters, global: undefined },
    { ...parameters, posterior: { ...parameters.posterior, lambda: 1.01 } },
    { ...parameters, frontal: { ...parameters.frontal, tau2: 100 } },
    ...[['U', -.1], ['beta', 22], ['beta', '11'], ['a', .5], ['tau1', 0], ['density', .3], ['cue', 49], ['cue', .5]].map(([key, value]) =>
      ({ ...parameters, global: { ...parameters.global, [key]: value } }))
  ]) assert.equal((await request(body)).status, 400);
  const html = await (await fetch(`${base}/models/potts_topology`)).text();
  assert.match(html, /Many-to-many Potts Network/);
  assert.equal((html.match(/<option value="T\d\d"/g) ?? []).length, 12);
  assert.match(html, /score-help-2/);
  console.log('PASS: fresh concurrent seeds and samples, paired 49-memory API, baseline and many-to-many choices, global controls, validation, third-model page.');
} finally { await server.close(); }
