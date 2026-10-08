import { json } from '@sveltejs/kit';
import type { RequestHandler } from './$types';
import { streamSimulation } from '$lib/server/simulation';

export const config = { runtime: 'nodejs22.x', maxDuration: 300 };

export const POST: RequestHandler = async ({ request }) => {
  let body;
  try { body = await request.json(); } catch { return json({ message: 'Invalid JSON.' }, { status: 400 }); }
  const { model, S, w, tau2, steps = 5000, posterior, frontal, lambda, topology = 'one-to-one' } = body ?? {};
  const number = (value: unknown, min: number, max: number) =>
    typeof value === 'number' && Number.isFinite(value) && value >= min && value <= max;
  const region = (value: { S?: number; w?: number; tau2?: number } | undefined, maxTau: number) =>
    value && Number.isInteger(value.S) && number(value.S, 3, 11) &&
    number(value.w, 0.6, 2) && number(value.tau2, 100, maxTau);
  if (!Number.isInteger(steps) || !number(steps, 1, 5000)) {
    return json({ message: 'Please choose between 1 and 5000 sweeps.' }, { status: 400 });
  }
  let args: number[];
  if (model === 'homo' && Number.isInteger(S) && number(S, 3, 11) && number(w, 0.6, 2) && number(tau2, 100, 800)) {
    args = [S, w, tau2, steps];
  } else if (model === 'hybrid' && region(posterior, 800) && region(frontal, 1600) &&
    number(lambda, 0, 1) && topology === 'one-to-one') {
    if (frontal.tau2 <= posterior.tau2) {
      return json({ message: 'Choose a slower frontal adaptation time than the posterior adaptation time.' }, { status: 400 });
    }
    args = [posterior.w, frontal.w, posterior.S, frontal.S, posterior.tau2, frontal.tau2, lambda, steps];
  } else {
    return json({ message: 'Please choose valid model parameters.' }, { status: 400 });
  }
  return new Response(streamSimulation(model, args, request.signal), {
    headers: { 'Content-Type': 'application/x-ndjson', 'Cache-Control': 'no-store', 'X-Accel-Buffering': 'no' }
  });
};
