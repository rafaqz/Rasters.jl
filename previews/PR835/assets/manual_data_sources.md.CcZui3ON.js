import{_ as a,c as s,o as t,aA as i}from"./chunks/framework.CeTlh2ns.js";const r="/Rasters.jl/previews/PR835/assets/qzuyden.Y4d2sml7.png",k=JSON.parse('{"title":"Data sources","description":"","frontmatter":{},"headers":[],"relativePath":"manual/data_sources.md","filePath":"manual/data_sources.md","lastUpdated":null}'),o={name:"manual/data_sources.md"};function n(d,e,l,c,h,p){return t(),s("div",null,e[0]||(e[0]=[i(`<h1 id="Data-sources" tabindex="-1">Data sources <a class="header-anchor" href="#Data-sources" aria-label="Permalink to &quot;Data sources {#Data-sources}&quot;">​</a></h1><p>Rasters.jl uses a number of backends to load raster data. <code>Raster</code>, <code>RasterStack</code> and <code>RasterSeries</code> will detect which backend to use for you, automatically.</p><h2 id="gdal" tabindex="-1">GDAL <a class="header-anchor" href="#gdal" aria-label="Permalink to &quot;GDAL&quot;">​</a></h2><p>All files GDAL can access, such as <code>.tiff</code> and <code>.asc</code> files, can be loaded, using <a href="https://github.com/yeesian/ArchGDAL.jl" target="_blank" rel="noreferrer">ArchGDAL.jl</a>. These are generally best loaded as <code>Raster(&quot;filename.tif&quot;)</code>, but can be loaded as <code>RasterStack(&quot;filename.tif&quot;; layersfrom=Band)</code>, taking layers from the <code>Band</code> dimension, which is also the default.</p><h2 id="netcdf" tabindex="-1">NetCDF <a class="header-anchor" href="#netcdf" aria-label="Permalink to &quot;NetCDF&quot;">​</a></h2><p>NetCDF <code>.nc</code> and some HDF5 <code>.h5</code> files cab be loaded using <a href="https://github.com/Alexander-Barth/NCDatasets.jl" target="_blank" rel="noreferrer">NCDatasets.jl</a>. Layers from files can be loaded as <code>Raster(&quot;filename.nc&quot;; name=:layername)</code>. Without <code>name</code> the first layer is used. <code>RasterStack(&quot;filename.nc&quot;)</code> will use all netcdf variables in the file that are not dimensions as layers.</p><p>NetCDF layers can have arbitrary dimensions. Known, common dimension names are converted to <code>X</code>, <code>Y</code> <code>Z</code>, and <code>Ti</code>, otherwise <code>Dim{:layername}</code> is used. Layers in the same file may also have different dimensions.</p><h2 id="zarr" tabindex="-1">Zarr <a class="header-anchor" href="#zarr" aria-label="Permalink to &quot;Zarr&quot;">​</a></h2><p>Zarr files can be loaded with the <a href="https://github.com/JuliaGeo/ZarrDatasets.jl" target="_blank" rel="noreferrer">ZarrDatasets.jl</a> backend. <code>Raster(filename; source=Zarrsource())</code> may be needed where the file type cant be detected from the filename. <code>write</code> does not yet work for Zarr but will in future.</p><h2 id="grib" tabindex="-1">GRIB <a class="header-anchor" href="#grib" aria-label="Permalink to &quot;GRIB&quot;">​</a></h2><p>GRIB files can be loaded with the <a href="https://github.com/JuliaGeo/GRIBDatasets.jl" target="_blank" rel="noreferrer">ZarrDatasets.jl</a>. <code>write</code> is not implemented for GRIB.</p><h2 id="grd" tabindex="-1">GRD <a class="header-anchor" href="#grd" aria-label="Permalink to &quot;GRD&quot;">​</a></h2><p>R GRD files can be loaded natively, using Julias <code>MMap</code> - which means they are very fast, but are not compressed. They are always 3 dimensional, and have <code>Y</code>, <code>X</code> and <a href="/Rasters.jl/previews/PR835/api#Rasters.Band"><code>Band</code></a> dimensions.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Rasters</span></span></code></pre></div><h2 id="Writing-file-formats-to-disk" tabindex="-1">Writing file formats to disk <a class="header-anchor" href="#Writing-file-formats-to-disk" aria-label="Permalink to &quot;Writing file formats to disk {#Writing-file-formats-to-disk}&quot;">​</a></h2><p>Files can be written to disk with ArchGDAL.jl and NCDatasets.jl backends using <code>write(&quot;filename.ext&quot;, raster)</code>. See the docs for <a href="/Rasters.jl/previews/PR835/api#Base.write-Tuple{AbstractString, AbstractRasterSeries}"><code>write</code></a>.</p><p>They can (with some caveats) be written to different formats than they were loaded in as, providing file-type conversion for spatial data. Some metadata may be lost in formats that store little metadata, or where metadata conversion has not been completely implemented.</p><h2 id="RasterDataSources.jl-integration" tabindex="-1">RasterDataSources.jl integration <a class="header-anchor" href="#RasterDataSources.jl-integration" aria-label="Permalink to &quot;RasterDataSources.jl integration {#RasterDataSources.jl-integration}&quot;">​</a></h2><p><a href="https://github.com/EcoJulia/RasterDataSources.jl" target="_blank" rel="noreferrer">RasterDataSources.jl</a> standardises the download of common raster data sources, with a focus on datasets used in ecology and the environmental sciences. RasterDataSources.jl is tightly integrated into Rasters.jl, so that datsets and keywords can be used directly to download and load data as a <code>Raster</code>, <code>RasterStack</code>, or <code>RasterSeries</code>.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Rasters, CairoMakie, Dates</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> RasterDataSources</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">A </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Raster</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(WorldClim{Climate}, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:tavg</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">; month</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">June)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">Makie</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">plot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(A)</span></span></code></pre></div><p><img src="`+r+'" alt="" width="600px" height="450px"></p><p>See the docs for <a href="/Rasters.jl/previews/PR835/api#Rasters.Raster"><code>Raster</code></a>, <a href="/Rasters.jl/previews/PR835/api#Rasters.RasterStack"><code>RasterStack</code></a> and <a href="/Rasters.jl/previews/PR835/api#Rasters.RasterSeries"><code>RasterSeries</code></a>, and the docs for <code>RasterDataSources.getraster</code> for syntax to specify various data sources.</p>',22)]))}const f=a(o,[["render",n]]);export{k as __pageData,f as default};
