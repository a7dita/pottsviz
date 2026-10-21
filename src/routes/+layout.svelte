<script lang="ts">
	import 'virtual:uno.css';
	import '@unocss/reset/tailwind.css';
</script>

<svelte:head>
	<title>PottsViz — Explore memory latching</title>
	<meta name="description" content="An interactive toy app for exploring how memories persist, fade, or give way to others in Potts associative neural networks." />
	<meta property="og:type" content="website" />
	<meta property="og:site_name" content="PottsViz" />
	<meta property="og:title" content="PottsViz — Explore memory latching" />
	<meta property="og:description" content="Change connections and dynamics, and watch how one memory leads to another in Potts associative neural networks." />
	<meta property="og:image" content="https://pottsviz.vercel.app/potts-homo-share.jpg" />
	<meta property="og:image:type" content="image/jpeg" />
	<meta property="og:image:width" content="1200" />
	<meta property="og:image:height" content="951" />
	<meta property="og:image:alt" content="Homogeneous Potts network: interconnected units with coloured states over a brain illustration." />
	<meta name="twitter:card" content="summary_large_image" />
	<meta name="twitter:title" content="PottsViz — Explore memory latching" />
	<meta name="twitter:description" content="Explore how memories persist, fade, or give way to others in an interactive Potts neural network." />
	<meta name="twitter:image" content="https://pottsviz.vercel.app/potts-homo-share.jpg" />
	<meta name="twitter:image:alt" content="Homogeneous Potts network: interconnected units with coloured states over a brain illustration." />
</svelte:head>

<div class="flex h-80vh text-gray-800 font-sans">
	<div class="m-auto">
		<slot />
	</div>
</div>
