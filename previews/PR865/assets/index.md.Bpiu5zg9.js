import{_ as a,c as i,ai as e,o as t}from"./chunks/framework.DFAlnjdJ.js";const g=JSON.parse(`{"title":"","description":"","frontmatter":{"layout":"home","hero":{"name":"Rasters.jl","text":"Manipulating spatial data","tagline":"a powerful package that simplifies the handling of rasterized spatial data in Julia.","image":{"src":"/logo.png","alt":"Rasters"},"actions":[{"theme":"brand","text":"Get Started","link":"/get_started"},{"theme":"alt","text":"View on Github","link":"https://github.com/rafaqz/Rasters.jl"},{"theme":"alt","text":"API Reference","link":"/api"}]},"features":[{"title":"💥💠 Core Functionality","details":"Defines common types and methods for reading, writing, and manipulating rasterized spatial data. <a class=\\"highlight-link\\">Rasters.jl</a> provides unified data handling through types like <a class=\\"highlight-link\\">Raster</a>, <a class=\\"highlight-link\\">RasterStack</a>, and <a class=\\"highlight-link\\">RasterSeries</a>, offering seamless abstraction regardless of storage backend.","link":"/manual/methods"},{"title":"⚙️ Data-Source Abstraction","details":"Rasters provides a standardized interface that enables various data source types to be used with consistent syntax. The data can include <a class=\\"highlight-link\\">GeoTIFF</a> or <a class=\\"highlight-link\\">NetCDF</a> files, in-memory Arrays, or <a class=\\"highlight-link\\">CuArrays</a> on the <a class=\\"highlight-link\\">GPU</a>—all of which will behave in the same way.","link":"/manual/data_sources"},{"title":"🗂️🌐 Data Formats","details":"A <a class=\\"highlight-link\\">RasterStack</a> can be backed by a <a class=\\"highlight-link\\">NetCDF</a> or <a class=\\"highlight-link\\">HDF5</a> file, a NamedTuple of Rasters holding <a class=\\"highlight-link\\">.tif</a> files, or all Rasters in memory. Users do not need to worry about the specifics of spatial file types."},{"title":"🌍🔍 Effortless Spatial Lookups","details":"<a class=\\"highlight-link\\">Projected</a> lookups with cylindrical projections can be indexed using other cylindrical projections by setting the <a class=\\"highlight-link\\">mappedcrs</a> keyword during construction. You don’t need to know the underlying projection, as the conversion is handled automatically. This means that <a class=\\"highlight-link\\">lat/lon EPSG(4326)</a> can be used seamlessly if needed."},{"title":"🧩⚡ Modular Extensions with Faster Loading","details":"From version 0.8, Rasters.jl uses Julia's package extension system. On Julia 1.9, packages can be placed in <a class=\\"highlight-link\\">extensions</a> to load only when needed, reducing loading time to under a second. However, each backend or feature <a class=\\"highlight-link\\">must be loaded manually</a>.","link":"/#Package-extensions"},{"title":"🐞📌 Bug reports","details":"Raster data is complex, and there are many opportunities for subtle or obvious bugs to appear. We need bug reports to reduce how often these issues occur over time. However, it's also crucial that issues are easy to reproduce; otherwise, it becomes impractical to fix them. Please follow <a href=\\"#Bugs\\" class=\\"highlight-link\\">these guidelines</a>!","link":"/#Bugs"}]},"headers":[],"relativePath":"index.md","filePath":"index.md","lastUpdated":null}`),l={name:"index.md"};function n(h,s,p,o,r,d){return t(),i("div",null,s[0]||(s[0]=[e(`<h2 id="How-to-Install-Rasters.jl?" tabindex="-1">How to Install Rasters.jl? <a class="header-anchor" href="#How-to-Install-Rasters.jl?" aria-label="Permalink to &quot;How to Install Rasters.jl? {#How-to-Install-Rasters.jl?}&quot;">​</a></h2><p>Since <code>Rasters.jl</code> is registered in the Julia General registry, you can simply run the following command in the Julia REPL:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">julia</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Pkg</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">julia</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Pkg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">add</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Rasters.jl&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># or</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">julia</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> ] </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># &#39;]&#39; should be pressed</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">pkg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> add Rasters</span></span></code></pre></div><p>If you want to use the latest unreleased version, you can run the following command:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">pkg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> add Rasters</span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;">#main</span></span></code></pre></div><h2 id="Package-extensions" tabindex="-1">Package extensions <a class="header-anchor" href="#Package-extensions" aria-label="Permalink to &quot;Package extensions {#Package-extensions}&quot;">​</a></h2><p>Before using a specific backend or feature, you must install and load its corresponding package,</p><div class="vp-code-group vp-adaptive-theme"><div class="tabs"><input type="radio" name="group-UCrAJ" id="tab-a8rLoZB" checked><label data-title=" gdal " for="tab-a8rLoZB"> gdal </label><input type="radio" name="group-UCrAJ" id="tab-yc8kB8r"><label data-title=" netcdf " for="tab-yc8kB8r"> netcdf </label><input type="radio" name="group-UCrAJ" id="tab-AjqG4_g"><label data-title=" grd " for="tab-AjqG4_g"> grd </label><input type="radio" name="group-UCrAJ" id="tab-WrwpbiQ"><label data-title=" smap " for="tab-WrwpbiQ"> smap </label><input type="radio" name="group-UCrAJ" id="tab-keH8cIc"><label data-title=" grib " for="tab-keH8cIc"> grib </label><input type="radio" name="group-UCrAJ" id="tab-A-eKOVL"><label data-title=" Makie plots " for="tab-A-eKOVL"> Makie plots </label><input type="radio" name="group-UCrAJ" id="tab-mmvQkgZ"><label data-title=" Data downloads " for="tab-mmvQkgZ"> Data downloads </label><input type="radio" name="group-UCrAJ" id="tab-BenDKYC"><label data-title=" Coordinate transforms " for="tab-BenDKYC"> Coordinate transforms </label></div><div class="blocks"><div class="language-julia vp-adaptive-theme active"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Pkg</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Pkg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">add</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;ArchGDAL&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Pkg</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Pkg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">add</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;NCDatasets&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># built-in</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Pkg</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Pkg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">add</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;HDF5&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Pkg</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Pkg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">add</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;GRIBDatasets&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Pkg</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Pkg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">add</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">([</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;GLMakie&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">])</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Pkg</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Pkg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">add</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">([</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;RasterDataSources&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]) </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># which includes Worldclim{Climate}</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Pkg</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Pkg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">add</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">([</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;CoordinateTransformations&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">]) </span><span style="--shiki-light:#6A737D;--shiki-dark:#6A737D;"># transformations for gdal rasters</span></span></code></pre></div></div></div><p>and as an example, to use the GDAL backend and download files, you will need to do:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Rasters, ArchGDAL, RasterDataSources</span></span></code></pre></div><h2 id="Bugs" tabindex="-1">🐞📌 Bugs <a class="header-anchor" href="#Bugs" aria-label="Permalink to &quot;🐞📌  Bugs {#Bugs}&quot;">​</a></h2><div class="info custom-block"><p class="custom-block-title">Bugs, errors and making issues for Rasters.jl</p><p>Because there are so many raster file types and variations of them, most of the time we need the <code>exact file</code> that caused your problem to know how to fix it, and be sure that we have actually fixed it when we are done. So fixing a Rasters.jl bug nearly always involves downloading some file and running some code that breaks with it (if you can trigger the bug without a file, that&#39;s great! but its not always possible).</p><p>To make an issue we can fix quickly (or at all) there are three key steps:</p><ol><li><p>Include the file in an accessible place on web <code>without authentication</code> or any other work on our part, so we can just get it and find your bug. You can put it on a file hosting platform (e.g. google drive, drop box, whatever you use) and share the url.</p></li><li><p>Add a <a href="https://discourse.julialang.org/t/please-read-make-it-easier-to-help-you/14757" target="_blank" rel="noreferrer"><code>minimum working example</code></a> to the issue template that first downloads the file, then runs the function that triggers the bug.</p></li><li><p>Paste the <code>complete stack trace</code> of the error it produces, right to the bottom, into the issue template. Then we can be sure we reproduced the same problem.</p></li></ol><p>Good issues are really appreciated, but they do take just a little extra effort with Rasters.jl because of this need for files.</p></div>`,12)]))}const c=a(l,[["render",n]]);export{g as __pageData,c as default};
