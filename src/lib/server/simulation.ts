import { spawn } from 'node:child_process';
import { resolve } from 'node:path';

/** One child process and one response per visitor; no shared files or jobs. */
export function streamSimulation(program: 'homo', args: number[], signal: AbortSignal) {
  const encoder = new TextEncoder();
  let stop = () => {};
  return new ReadableStream<Uint8Array>({
    start(controller) {
      const child = spawn(resolve('.simulators', program), [
        ...args.slice(0, 3).map(String), '--stream', String(args[3])
      ], { stdio: ['ignore', 'pipe', 'pipe'] });
      let ended = false;
      let pending = '';
      let diagnostics = '';
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
      }, 55000);
      signal.addEventListener('abort', stop, { once: true });
      send({ type: 'started', steps: args[3], memories: 100 });
      child.stdout.setEncoding('utf8');
      child.stdout.on('data', (chunk: string) => {
        pending += chunk;
        const lines = pending.split('\n');
        pending = lines.pop() ?? '';
        for (const line of lines) {
          const cells = line.trim().split(/\s+/).map(Number);
          if (cells.length === 101 && cells.every(Number.isFinite)) {
            send({ type: 'sample', time: cells[0], posterior: cells.slice(1) });
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
        } else finish({ type: 'done' });
      });
      if (signal.aborted) stop();
    },
    cancel() { stop(); }
  });
}
