import { spawn } from 'node:child_process';
import { resolve } from 'node:path';

/** One child process and one response per visitor; no shared files or jobs. */
export function streamSimulation(program: 'homo' | 'hybrid' | 'topology', args: (number | string)[], signal: AbortSignal, metadata: { memories?: number; topology?: string; patternSeed?: number; runtimeSeed?: number } = {}) {
  const encoder = new TextEncoder();
  let stop = () => {};
  return new ReadableStream<Uint8Array>({
    start(controller) {
      const steps = args[args.length - 1];
      const child = spawn(resolve('.simulators', program), program === 'homo' ? [
        ...args.slice(0, 3).map(String), '--stream', String(steps)
      ] : args.map(String), { stdio: ['ignore', 'pipe', 'pipe'] });
      let ended = false;
      let pending = '';
      let diagnostics = '';
      let posterior: number[] | undefined;
      let samples = 0;
      const startedAt = Date.now();
      const send = (event: object) => {
        if (!ended) controller.enqueue(encoder.encode(JSON.stringify(event) + '\n'));
      };
      const finish = (event?: object) => {
        if (ended) return;
        if (event) send(event);
        ended = true;
        clearTimeout(deadline);
        signal.removeEventListener('abort', stop);
        controller.close();
      };
      stop = () => {
        if (ended) return;
        console.info('Simulation cancelled', { program, samples, elapsedMs: Date.now() - startedAt });
        child.kill('SIGKILL');
        if (!ended) {
          ended = true;
          clearTimeout(deadline);
          signal.removeEventListener('abort', stop);
        }
      };
      const deadline = setTimeout(() => {
        child.kill('SIGKILL');
        finish({ type: 'error', message: 'Simulation time limit reached. Please try a shorter run.' });
      }, program === 'homo' ? 55000 : 270000);
      signal.addEventListener('abort', stop, { once: true });
      const memories = metadata.memories ?? 100;
      send({ type: 'started', steps, memories, topology: program === 'hybrid' ? 'one-to-one' : undefined, ...metadata });
      child.stdout.setEncoding('utf8');
      child.stdout.on('data', (chunk: string) => {
        pending += chunk;
        const lines = pending.split('\n');
        pending = lines.pop() ?? '';
        for (const line of lines) {
          const cells = line.trim().split(/\s+/).map(Number);
          if (cells.length === memories + 1 && cells.every(Number.isFinite)) {
            if (program === 'homo') { samples++; send({ type: 'sample', time: cells[0], posterior: cells.slice(1) }); }
            else if (!posterior) posterior = cells;
            else {
              if (posterior[0] !== cells[0]) {
                child.kill('SIGKILL');
                finish({ type: 'error', message: 'The network snapshots could not be synchronized.' });
                return;
              }
              send({ type: 'sample', time: cells[0], posterior: posterior.slice(1), frontal: cells.slice(1) });
              samples++;
              posterior = undefined;
            }
          }
        }
      });
      child.stderr.on('data', (chunk) => { diagnostics = (diagnostics + String(chunk)).slice(-2000); });
      child.on('error', (error) => {
        console.error('Simulator failed to launch', error);
        finish({ type: 'error', message: 'The simulation could not start. Please try again.' });
      });
      child.on('close', (code) => {
        if (ended) return;
        if (code !== 0) {
          console.error('Simulator exited', { program, code, diagnostics });
          finish({ type: 'error', message: 'The simulation failed. Please try again.' });
        } else if (posterior) finish({ type: 'error', message: 'The final network snapshot was incomplete.' });
        else {
          console.info('Simulation finished', { program, samples, elapsedMs: Date.now() - startedAt });
          finish(samples ? { type: 'done' } : { type: 'error', message: 'The simulation produced no valid snapshots.' });
        }
      });
      if (signal.aborted) stop();
    },
    cancel() { stop(); }
  });
}
