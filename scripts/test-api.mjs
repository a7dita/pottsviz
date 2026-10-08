import assert from 'node:assert/strict';
import { createServer } from 'vite';
import { readFileSync } from 'node:fs';
import { execFileSync } from 'node:child_process';

const reference = JSON.parse(readFileSync('tests/ryom-reference.json', 'utf8'));
for (const fixture of reference.cases) {
  const output = execFileSync('.simulators/hybrid', fixture.args, { encoding: 'utf8', stdio: ['ignore', 'pipe', 'ignore'] });
  const rows = output.trim().split('\n').map(line => line.trim().split(/\s+/).map(Number));
  assert.deepEqual(rows, fixture.rows, 'Original Ryom dynamics and one-to-one coupling are preserved');
}
console.log('PASS: three numerical regression fixtures from the original Ryom implementation.');

const server = await createServer({ server: { host: '127.0.0.1', port: 0 } });
await server.listen();
const address = server.httpServer.address();
const base = `http://127.0.0.1:${address.port}`;
const parameters = { model: 'homo', S: 7, w: 1.2, tau2: 200, steps: 30 };
async function request(body, path = '/api/simulate') {
  return fetch(base + path, {
    method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body)
  });
}
async function simulate(body) {
  const response = await request(body);
  assert.equal(response.status, 200);
  assert.match(response.headers.get('content-type'), /ndjson/);
  const events = (await response.text()).trim().split('\n').map(line => JSON.parse(line));
  assert.equal(events[0].type, 'started');
  assert.equal(events.at(-1).type, 'done');
  const samples = events.filter(event => event.type === 'sample');
  assert.equal(samples.length, 7);
  assert.equal(samples.at(-1).time, 30);
  assert(samples.every(sample => sample.posterior.length === 100 && sample.posterior.every(Number.isFinite)));
  return samples;
}
try {
  const [first, second] = await Promise.all([simulate(parameters), simulate(parameters)]);
  assert.deepEqual(first, second, 'Independent simultaneous visitors get the same seeded result');
  await simulate({ ...parameters, S: 3, w: 0.6, tau2: 800 });
  const hybrid = {
    model: 'hybrid', topology: 'one-to-one', lambda: 0.9, steps: 30,
    posterior: { S: 3, w: 1.1, tau2: 100 }, frontal: { S: 7, w: 1.1, tau2: 400 }
  };
  const response = await request(hybrid);
  assert.equal(response.status, 200);
  const events = (await response.text()).trim().split('\n').map(line => JSON.parse(line));
  assert.equal(events.at(-1).type, 'done');
  const pairs = events.filter(event => event.type === 'sample');
  assert.equal(pairs.length, 4);
  assert.equal(pairs.at(-1).time, 30);
  assert(pairs.every(sample => sample.posterior.length === 100 && sample.frontal.length === 100 &&
    [...sample.posterior, ...sample.frontal].every(Number.isFinite)));
  for (const invalid of [
    { ...hybrid, topology: 'many-to-many' }, { ...hybrid, lambda: 1.1 },
    { ...hybrid, frontal: { ...hybrid.frontal, tau2: 100 } },
    { ...hybrid, posterior: { ...hybrid.posterior, S: 12 } }
  ]) assert.equal((await request(invalid)).status, 400);
  for (const invalid of [{ ...parameters, S: '7' }, { ...parameters, w: null },
    { ...parameters, steps: 5001 }, { command: 'shell commands are not model parameters' }]) {
    assert.equal((await request(invalid)).status, 400);
  }
  for (const path of ['/api', '/api2']) assert.equal((await request({}, path)).status, 410);
  console.log('PASS: homogeneous and paired fronto-posterior live output, independent concurrent runs, parameter bounds, retired endpoints.');
} finally { await server.close(); }
