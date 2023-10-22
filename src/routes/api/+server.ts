import { json } from '@sveltejs/kit';
import type { RequestHandler } from './$types';
import { exec } from 'node:child_process';
import fs from 'node:fs/promises';
import util from 'node:util';


const readOutput = async (textfileName: string) => {
	try {
		const output = await fs.readFile(textfileName, 'utf8');
		return output;
	} catch (err) {
		return err
	}
}

export const POST: RequestHandler = async ({ request }) => {

	const { textfileName } = await request.json();
	const output = await readOutput(textfileName);
	return new Response(output);

};
