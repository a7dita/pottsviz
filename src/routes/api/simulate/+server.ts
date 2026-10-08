import { json } from '@sveltejs/kit';
import type { RequestHandler } from './$types';
import { streamSimulation } from '$lib/server/simulation';

export const config = { runtime: 'nodejs22.x', maxDuration: 60 };

export const POST: RequestHandler = async ({ request }) => {
  let body;
  try { body = await request.json(); } catch { return json({ message: 'Invalid JSON.' }, { status: 400 }); }
  const { model, S, w, tau2, steps = 5000 } = body ?? {};
  if (model !== 'homo' || !Number.isInteger(S) || S < 3 || S > 11 ||
      typeof w !== 'number' || !Number.isFinite(w) || w < 0.6 || w > 2 ||
      typeof tau2 !== 'number' || !Number.isFinite(tau2) || tau2 < 100 || tau2 > 800 ||
      !Number.isInteger(steps) || steps < 1 || steps > 5000) {
    return json({ message: 'Please choose valid model parameters.' }, { status: 400 });
  }
  return new Response(streamSimulation('homo', [S, w, tau2, steps], request.signal), {
    headers: { 'Content-Type': 'application/x-ndjson', 'Cache-Control': 'no-store', 'X-Accel-Buffering': 'no' }
  });
};
