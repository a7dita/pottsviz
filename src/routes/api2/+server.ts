import { json } from '@sveltejs/kit';
export const POST = () => json({ message: 'Please use /api/simulate with model parameters.' }, { status: 410 });
