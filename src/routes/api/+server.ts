import { json } from '@sveltejs/kit';
import type { RequestHandler } from './$types';
import { exec } from 'node:child_process';
import fs from 'node:fs/promises';
import util from 'node:util';

// this helper function reads the textfile
const readOutput = async (textfileName: string) => {
	try {
		const output = await fs.readFile(textfileName, 'utf8');
		return output;
	} catch (err) {
		return err
	}
}

// upon getting a POST request with the textfile name, this api reads the content of the file
// and send it (as a text string) to the frontend.
export const POST: RequestHandler = async ({ request }) => {

	const { textfileName } = await request.json();
	const output = await readOutput(textfileName);
	return new Response(output);

};
