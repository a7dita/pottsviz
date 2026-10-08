export interface Sample { time: number; posterior: number[]; frontal?: number[] }

export interface RunInfo { memories: number; topology?: string; patternSeed?: number; runtimeSeed?: number; posteriorCue?: number | 'random'; frontalCue?: number | 'random' }

export async function runSimulation(
  parameters: object, signal: AbortSignal, onSample: (sample: Sample) => void, onStarted?: (info: RunInfo) => void
) {
  const response = await fetch('/api/simulate', {
    method: 'POST', headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(parameters), signal
  });
  if (!response.ok) {
    const error = await response.json().catch(() => ({}));
    throw new Error(error.message ?? `Simulation request failed (${response.status}).`);
  }
  if (!response.body) throw new Error('This browser could not read the live simulation.');
  const reader = response.body.getReader();
  const decoder = new TextDecoder();
  let pending = '';
  let completed = false;
  try {
    while (true) {
      const { done, value } = await reader.read();
      pending += decoder.decode(value, { stream: !done });
      const lines = pending.split('\n');
      pending = lines.pop() ?? '';
      for (const line of lines) {
        if (!line.trim()) continue;
        const event = JSON.parse(line);
        if (event.type === 'started') onStarted?.(event);
        if (event.type === 'sample') onSample(event);
        if (event.type === 'error') throw new Error(event.message);
        if (event.type === 'done') completed = true;
      }
      if (done) break;
    }
    if (!completed) throw new Error('The simulation connection ended early. Please try again.');
  } finally {
    await reader.cancel().catch(() => {});
    reader.releaseLock();
  }
}
