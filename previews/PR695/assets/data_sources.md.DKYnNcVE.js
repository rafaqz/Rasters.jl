import{_ as a,c as s,a5 as t,o as i}from"./chunks/framework.JYPsFMF-.js";const o="/Rasters.jl/previews/PR695/assets/vitlztl.CBMwZ7qg.png",m=JSON.parse('{"title":"Data sources","description":"","frontmatter":{},"headers":[],"relativePath":"data_sources.md","filePath":"data_sources.md","lastUpdated":null}'),r={name:"data_sources.md"};function n(d,e,l,c,h,p){return i(),s("div",null,e[0]||(e[0]=[t(`<h1 id="Data-sources" tabindex="-1">Data sources <a class="header-anchor" href="#Data-sources" aria-label="Permalink to &quot;Data sources {#Data-sources}&quot;">​</a></h1><p>Rasters.jl uses a number of backends to load raster data. <code>Raster</code>, <code>RasterStack</code> and <code>RasterSeries</code> will detect which backend to use for you, automatically.</p><h2 id="grd" tabindex="-1">GRD <a class="header-anchor" href="#grd" aria-label="Permalink to &quot;GRD&quot;">​</a></h2><p>R GRD files can be loaded natively, using Julias <code>MMap</code> - which means they are very fast, but are not compressed. They are always 3 dimensional, and have <code>Y</code>, <code>X</code> and <a href="/Rasters.jl/previews/PR695/api#Rasters.Band"><code>Band</code></a> dimensions.</p><h2 id="netcdf" tabindex="-1">NetCDF <a class="header-anchor" href="#netcdf" aria-label="Permalink to &quot;NetCDF&quot;">​</a></h2><p>NetCDF <code>.nc</code> files are loaded using <a href="https://github.com/Alexander-Barth/NCDatasets.jl" target="_blank" rel="noreferrer">NCDatasets.jl</a>. Layers from files can be loaded as <code>Raster(&quot;filename.nc&quot;; name=:layername)</code>. Without <code>name</code> the first layer is used. <code>RasterStack(&quot;filename.nc&quot;)</code> will use all netcdf variables in the file that are not dimensions as layers.</p><p>NetCDF layers can have arbitrary dimensions. Known, common dimension names are converted to <code>X</code>, <code>Y</code> <code>Z</code>, and <code>Ti</code>, otherwise <code>Dim{:layername}</code> is used. Layers in the same file may also have different dimensions.</p><p>NetCDF files still have issues loading directly from disk for some operations. Using <code>read(ncstack)</code> may help.</p><h2 id="gdal" tabindex="-1">GDAL <a class="header-anchor" href="#gdal" aria-label="Permalink to &quot;GDAL&quot;">​</a></h2><p>All files GDAL can access, such as <code>.tiff</code> and <code>.asc</code> files, can be loaded, using <a href="https://github.com/yeesian/ArchGDAL.jl/issues" target="_blank" rel="noreferrer">ArchGDAL.jl</a>. These are generally best loaded as <code>Raster(&quot;filename.tif&quot;)</code>, but can be loaded as <code>RasterStack(&quot;filename.tif&quot;; layersfrom=Band)</code>, taking layers from the <code>Band</code> dimension, which is also the default.</p><h2 id="smap" tabindex="-1">SMAP <a class="header-anchor" href="#smap" aria-label="Permalink to &quot;SMAP&quot;">​</a></h2><p>The <a href="https://smap.jpl.nasa.gov/" target="_blank" rel="noreferrer">Soil Moisture Active-Passive</a> dataset provides global layers of soil moisture, temperature and other related data, in a custom HDF5 format. Layers are always 2 dimensional, with <code>Y</code> and <code>X</code> dimensions.</p><p>These can be loaded as multi-layered <code>RasterStack(&quot;filename.h5&quot;)</code>. Individual layers can be loaded as <code>Raster(&quot;filename.h5&quot;; name=:layername)</code>, without <code>name</code> the first layer is used.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Rasters</span></span></code></pre></div><div class="warning custom-block"><p class="custom-block-title">Missing docstring.</p><p>Missing docstring for <code>smapseries</code>. Check Documenter&#39;s build log for details.</p></div><h2 id="Writing-file-formats-to-disk" tabindex="-1">Writing file formats to disk <a class="header-anchor" href="#Writing-file-formats-to-disk" aria-label="Permalink to &quot;Writing file formats to disk {#Writing-file-formats-to-disk}&quot;">​</a></h2><p>Files can be written to disk in all formats other than SMAP HDF5 using <code>write(&quot;filename.ext&quot;, A)</code>. See the docs for <a href="/Rasters.jl/previews/PR695/api#Base.write-Tuple{AbstractString, AbstractRasterSeries}"><code>write</code></a>. They can (with some caveats) be written to different formats than they were loaded in as, providing file-type conversion for spatial data.</p><p>Some metadata may be lost in formats that store little metadata, or where metadata conversion has not been completely implemented.</p><h2 id="RasterDataSources.jl-integration" tabindex="-1">RasterDataSources.jl integration <a class="header-anchor" href="#RasterDataSources.jl-integration" aria-label="Permalink to &quot;RasterDataSources.jl integration {#RasterDataSources.jl-integration}&quot;">​</a></h2><p><a href="https://github.com/EcoJulia/RasterDataSources.jl" target="_blank" rel="noreferrer">RasterDataSources.jl</a> standardises the download of common raster data sources, with a focus on datasets used in ecology and the environmental sciences. RasterDataSources.jl is tightly integrated into Rasters.jl, so that datsets and keywords can be used directly to download and load data as a <code>Raster</code>, <code>RasterStack</code>, or <code>RasterSeries</code>.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Rasters, CairoMakie, Dates</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> RasterDataSources</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">A </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Raster</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(WorldClim{Climate}, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:tavg</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">; month</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">June)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Makie</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(A)</span></span></code></pre></div><p><img src="`+o+'" alt=""></p><p>See the docs for <a href="/Rasters.jl/previews/PR695/api#Rasters.Raster"><code>Raster</code></a>, <a href="/Rasters.jl/previews/PR695/api#Rasters.RasterStack"><code>RasterStack</code></a> and <a href="/Rasters.jl/previews/PR695/api#Rasters.RasterSeries"><code>RasterSeries</code></a>, and the docs for <code>RasterDataSources.getraster</code> for syntax to specify various data sources.</p>',23)]))}const k=a(r,[["render",n]]);export{m as __pageData,k as default};