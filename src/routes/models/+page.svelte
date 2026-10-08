<div class="space-y-4 flex flex-col p-10">
	<p class="text-purple text-4xl">The models you can try</p>
	<div class="inline-flex items-center justify-center w-full">
		<hr class="w-full h-px bg-gray-200 border-0 dark:bg-gray-700" />
	</div>

	<a href="/models/potts_homo">
		<div
			class="space-y-8 animate-pulse p-2 md:space-y-0 md:space-x-8 rtl:space-x-reverse md:flex md:items-center border-white border-6 hover:border-blue-200 hover-border-6"
		>
			<div
				class="flex items-center justify-center w-full bg-gray-100 rounded sm:w-96 dark:bg-gray-700"
			>
				<img class="p-2" src="/potts_homo.png" alt="Homogeneous Potts network" />
			</div>
			<div class="w-full p-6">
				<div class="mb-4 text-gray-600 text-2xl">
					<p>Homogenous Potts Network</p>
				</div>
				<div class="h-2 bg-gray-200 rounded-full dark:bg-gray-700 max-w-[480px] mb-2.5"></div>
				<div class="h-2 bg-gray-200 rounded-full dark:bg-gray-700 mb-2.5"></div>
				<div class="h-2 bg-gray-200 rounded-full dark:bg-gray-700 max-w-[440px] mb-2.5"></div>
				<div class="h-2 bg-gray-200 rounded-full dark:bg-gray-700 max-w-[460px] mb-2.5"></div>
				<div class="h-2 bg-gray-200 rounded-full dark:bg-gray-700 max-w-[260px]"></div>
			</div>
		</div>
	</a>

	<a href="/models/potts_hybrid">
		<div
			class="space-y-8 animate-pulse p-2 md:space-y-0 md:space-x-8 rtl:space-x-reverse md:flex md:items-center border-white border-6 hover:border-blue-200 hover-border-6"
		>
			<div
				class="flex items-center justify-center w-full bg-gray-100 rounded sm:w-96 dark:bg-gray-700"
			>
				<img class="p-2" src="/potts_hybrid.png" alt="Fronto-posterior Potts network" />
			</div>
			<div class="w-full p-6">
				<div class="mb-4 text-gray-600 text-2xl">
					<p>Hybrid Potts Network</p>
				</div>
				<div class="h-2 bg-gray-200 rounded-full dark:bg-gray-700 max-w-[480px] mb-2.5"></div>
				<div class="h-2 bg-gray-200 rounded-full dark:bg-gray-700 mb-2.5"></div>
				<div class="h-2 bg-gray-200 rounded-full dark:bg-gray-700 max-w-[440px] mb-2.5"></div>
				<div class="h-2 bg-gray-200 rounded-full dark:bg-gray-700 max-w-[460px] mb-2.5"></div>
				<div class="h-2 bg-gray-200 rounded-full dark:bg-gray-700 max-w-[260px]"></div>
			</div>
		</div>
	</a>
	<a href="/models/potts_topology" class="flex flex-wrap items-center gap-8 border-2 border-purple-100 rounded p-6 hover:bg-purple-50">
        <svg viewBox="0 0 240 160" role="img" aria-label="Many-to-many memory associations" style="width:240px;max-width:100%">
            <rect width="240" height="160" rx="12" fill="#f5f3ff" />
            {#each [0,1,2,3,4,5,6] as i}
                {#each [0,1,2,3,4,5,6] as j}
                    <line x1="40" y1={20+i*20} x2="200" y2={20+j*20} stroke="#c4b5fd" stroke-width="1" />
                {/each}
                <circle cx="40" cy={20+i*20} r="5" fill="#6d28d9" />
                <circle cx="200" cy={20+i*20} r="5" fill="#0891b2" />
            {/each}
        </svg>
        <div><h2 class="text-2xl text-gray-600">Many-to-many Potts Network</h2>
            <p class="mt-3">49 memories per region · 12 frozen association topologies</p>
            <p class="mt-2">Explore redundancy, modularity and mixing with global and regional controls.</p>
        </div>
    </a>
</div>
