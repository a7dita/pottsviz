import { sveltekit } from '@sveltejs/kit/vite';
import UnoCSS from 'unocss/vite';
import { presetWebFonts } from 'unocss';
import { presetUno } from 'unocss';
import { presetIcons } from 'unocss';
import { defineConfig } from 'vite';

export default defineConfig({
	plugins: [
		sveltekit(),
		UnoCSS({
			safelist: ['grid-cols-4', 'grid-cols-2'],
			presets: [
				presetUno(),
				presetIcons(),
				presetWebFonts({
					provider: 'bunny',
					fonts: {
						sans: 'Inter',
						mono: 'JetBrains Mono',
						barrio: 'Barrio'
					}
				})
			]
		})
	]
});
