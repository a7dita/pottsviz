<svelte:head>
	<title>Pottsviz</title>
	<meta
		name="description"
		content="Interactive Potts associative-network visualizations for homogeneous and fronto-posterior memory-latching models."
	/>
	<meta property="og:type" content="website" />
	<meta property="og:url" content="https://pottsviz.vercel.app/" />
	<meta property="og:site_name" content="Pottsviz" />
	<meta property="og:title" content="Pottsviz" />
	<meta
		property="og:description"
		content="Play with homogeneous and fronto-posterior Potts associative-network models in your browser."
	/>
	<meta property="og:image" content="https://pottsviz.vercel.app/potts_homo.png" />
	<meta property="og:image:secure_url" content="https://pottsviz.vercel.app/potts_homo.png" />
	<meta property="og:image:type" content="image/png" />
	<meta property="og:image:width" content="2502" />
	<meta property="og:image:height" content="1982" />
	<meta
		property="og:image:alt"
		content="A homogeneous Potts associative network over a brain illustration."
	/>
	<meta name="twitter:card" content="summary_large_image" />
	<meta name="twitter:title" content="Pottsviz" />
	<meta
		name="twitter:description"
		content="Play with homogeneous and fronto-posterior Potts associative-network models in your browser."
	/>
	<meta name="twitter:image" content="https://pottsviz.vercel.app/potts_homo.png" />
</svelte:head>

<script lang="ts">
	import 'virtual:uno.css';
	import '@unocss/reset/tailwind.css';
</script>

<div class="flex h-80vh text-gray-800 font-sans">
	<div class="m-auto">
		<slot />
	</div>
</div>
