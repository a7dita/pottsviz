import { json } from '@sveltejs/kit';
import type { RequestHandler } from './$types';
import { exec } from 'node:child_process';
import fs from 'node:fs/promises';
import util from 'node:util';


const readOutput = async () => {
	try {
		const output = await fs.readFile('data/mall_cue0_seed1_d', 'utf8');
		return output;
	} catch (err) {
		return err
	}
}

export const GET: RequestHandler = async ({ request }) => {
	const output = await readOutput();
	return new Response(output);

};
