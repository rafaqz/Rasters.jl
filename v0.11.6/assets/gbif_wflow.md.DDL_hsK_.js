import{_ as s,c as a,o as i,a7 as n}from"./chunks/framework.CyWcZaab.js";const p="/Rasters.jl/v0.11.6/assets/eclwhob.CVwZE_Z9.png",y=JSON.parse('{"title":"","description":"","frontmatter":{},"headers":[],"relativePath":"gbif_wflow.md","filePath":"gbif_wflow.md","lastUpdated":null}'),t={name:"gbif_wflow.md"},l=n(`<p>Load occurrences for the Mountain Pygmy Possum using GBIF.jl</p><h2 id="Load-GBIF" tabindex="-1">Load GBIF <a class="header-anchor" href="#Load-GBIF" aria-label="Permalink to &quot;Load GBIF {#Load-GBIF}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Rasters, GBIF2</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> RasterDataSources</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">const</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> RS </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Rasters</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Rasters</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">records </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> GBIF2</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">occurrence_search</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;Burramys parvus&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">; limit</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">300</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>300-element GBIF2.Table{GBIF2.Occurrence, JSON3.Array{JSON3.Object, Vector{UInt8}, SubArray{UInt64, 1, Vector{UInt64}, Tuple{UnitRange{Int64}}, true}}}┌──────────────────────────┬─────────┬─────────┬─────────┬──────────┬───────────</span></span>
<span class="line"><span>│                 geometry │    year │   month │     day │  kingdom │   phylum ⋯</span></span>
<span class="line"><span>│ Tuple{Float64, Float64}? │  Int64? │  Int64? │  Int64? │  String? │  String? ⋯</span></span>
<span class="line"><span>├──────────────────────────┼─────────┼─────────┼─────────┼──────────┼───────────</span></span>
<span class="line"><span>│                  missing │    2021 │       1 │       6 │ Animalia │ Chordata ⋯</span></span>
<span class="line"><span>│      (148.333, -36.4333) │    2011 │      11 │      21 │ Animalia │ Chordata ⋯</span></span>
<span class="line"><span>│      (148.396, -36.3818) │    2016 │      11 │      15 │ Animalia │ Chordata ⋯</span></span>
<span class="line"><span>│                  missing │ missing │ missing │ missing │ Animalia │ Chordata ⋯</span></span>
<span class="line"><span>│                  missing │ missing │ missing │ missing │ Animalia │ Chordata ⋯</span></span>
<span class="line"><span>│      (148.391, -36.3036) │    2015 │      11 │      15 │ Animalia │ Chordata ⋯</span></span>
<span class="line"><span>│                  missing │ missing │ missing │ missing │ Animalia │ Chordata ⋯</span></span>
<span class="line"><span>│                  missing │ missing │ missing │ missing │ Animalia │ Chordata ⋯</span></span>
<span class="line"><span>│      (147.096, -36.9357) │    2020 │       2 │      10 │ Animalia │ Chordata ⋯</span></span>
<span class="line"><span>│      (148.329, -36.4317) │    2016 │       1 │       3 │ Animalia │ Chordata ⋯</span></span>
<span class="line"><span>│                  missing │ missing │ missing │ missing │ Animalia │ Chordata ⋯</span></span>
<span class="line"><span>│                  missing │ missing │ missing │ missing │ Animalia │ Chordata ⋯</span></span>
<span class="line"><span>│                  missing │ missing │ missing │ missing │ Animalia │ Chordata ⋯</span></span>
<span class="line"><span>│      (148.347, -36.5047) │    2012 │      11 │      22 │ Animalia │ Chordata ⋯</span></span>
<span class="line"><span>│      (148.236, -36.5249) │    2012 │      11 │      23 │ Animalia │ Chordata ⋯</span></span>
<span class="line"><span>│            ⋮             │    ⋮    │    ⋮    │    ⋮    │    ⋮     │    ⋮     ⋱</span></span>
<span class="line"><span>└──────────────────────────┴─────────┴─────────┴─────────┴──────────┴───────────</span></span>
<span class="line"><span>                                                 78 columns and 285 rows omitted</span></span></code></pre></div><h2 id="Extract-coordinates" tabindex="-1">Extract coordinates <a class="header-anchor" href="#Extract-coordinates" aria-label="Permalink to &quot;Extract coordinates {#Extract-coordinates}&quot;">​</a></h2><p>Extract the longitude/latitude value to a <code>Vector</code> of points (a <code>Tuple</code> counts as a <code>(x, y)</code> point in GeoInterface.jl):</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">coords </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> [(r</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">decimalLongitude, r</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">decimalLatitude) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> r </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> records </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">if</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> !</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">ismissing</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(r</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">decimalLatitude)]</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>256-element Vector{Tuple{Float64, Float64}}:</span></span>
<span class="line"><span> (148.332969, -36.433349)</span></span>
<span class="line"><span> (148.396453, -36.381847)</span></span>
<span class="line"><span> (148.391097, -36.30362)</span></span>
<span class="line"><span> (147.096394, -36.935687)</span></span>
<span class="line"><span> (148.328896, -36.431684)</span></span>
<span class="line"><span> (148.347186, -36.504673)</span></span>
<span class="line"><span> (148.235596, -36.524924)</span></span>
<span class="line"><span> (148.240881, -36.400058)</span></span>
<span class="line"><span> (148.4167, -36.35)</span></span>
<span class="line"><span> (148.3167, -36.4167)</span></span>
<span class="line"><span> ⋮</span></span>
<span class="line"><span> (148.391682, -36.373215)</span></span>
<span class="line"><span> (148.394749, -36.284565)</span></span>
<span class="line"><span> (148.333783, -36.432552)</span></span>
<span class="line"><span> (148.333783, -36.432552)</span></span>
<span class="line"><span> (148.377673, -36.418261)</span></span>
<span class="line"><span> (148.328025, -36.437709)</span></span>
<span class="line"><span> (148.398438, -36.382602)</span></span>
<span class="line"><span> (148.310064, -36.448374)</span></span>
<span class="line"><span> (148.259489, -36.490719)</span></span></code></pre></div><h2 id="Get-layer-/-Band" tabindex="-1">Get layer / Band <a class="header-anchor" href="#Get-layer-/-Band" aria-label="Permalink to &quot;Get layer / Band {#Get-layer-/-Band}&quot;">​</a></h2><p>Get BioClim layers and subset to south-east Australia</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">A </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> RasterStack</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(WorldClim{BioClim}, (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">3</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">7</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">12</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">se_aus </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> A[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">X</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">138</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> ..</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 155</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">), </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Y</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">40</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> ..</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> -</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">25</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">), RS</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">Band</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)]</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>╭────────────────────╮</span></span>
<span class="line"><span>│ 102×89 RasterStack │</span></span>
<span class="line"><span>├────────────────────┴─────────────────────────────────────────────────── dims ┐</span></span>
<span class="line"><span>  ↓ X Projected{Float64} LinRange{Float64}(138.0, 154.83333333333331, 102) ForwardOrdered Regular Intervals{Start},</span></span>
<span class="line"><span>  → Y Projected{Float64} LinRange{Float64}(-25.16666666666667, -39.83333333333333, 89) ReverseOrdered Regular Intervals{Start}</span></span>
<span class="line"><span>├────────────────────────────────────────────────────────────────────── layers ┤</span></span>
<span class="line"><span>  :bio1  eltype: Float32 dims: X, Y size: 102×89</span></span>
<span class="line"><span>  :bio3  eltype: Float32 dims: X, Y size: 102×89</span></span>
<span class="line"><span>  :bio7  eltype: Float32 dims: X, Y size: 102×89</span></span>
<span class="line"><span>  :bio12 eltype: Float32 dims: X, Y size: 102×89</span></span>
<span class="line"><span>├────────────────────────────────────────────────────────────────────── raster ┤</span></span>
<span class="line"><span>  extent: Extent(X = (138.0, 154.99999999999997), Y = (-39.83333333333333, -25.000000000000004))</span></span>
<span class="line"><span>  missingval: -3.4f38</span></span>
<span class="line"><span>  crs: GEOGCS[&quot;WGS 84&quot;,DATUM[&quot;WGS_1984&quot;,SPHEROID[&quot;WGS 84&quot;,6378137,298.257223563,AUTHORITY[&quot;EPSG&quot;,&quot;7030&quot;]],AUTHORITY[&quot;EPSG&quot;,&quot;6326&quot;]],PRIMEM[&quot;Greenwich&quot;,0,AUTHORITY[&quot;EPSG&quot;,&quot;8901&quot;]],UNIT[&quot;degree&quot;,0.0174532925199433,AUTHORITY[&quot;EPSG&quot;,&quot;9122&quot;]],AXIS[&quot;Latitude&quot;,NORTH],AXIS[&quot;Longitude&quot;,EAST],AUTHORITY[&quot;EPSG&quot;,&quot;4326&quot;]]</span></span>
<span class="line"><span>└──────────────────────────────────────────────────────────────────────────────┘</span></span></code></pre></div><p>Plot BioClim predictors and scatter occurrence points on all subplots</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Plots</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">p </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> plot</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(se_aus);</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">kw </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (legend</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:none</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, opacity</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.5</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, markershape</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:cross</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, markercolor</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:black</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">foreach</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(i </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">-&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> scatter!</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(p, coords; subplot</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">i, kw</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">...</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">), </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">4</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">p</span></span></code></pre></div><p><img src="`+p+`" alt=""></p><p>Then extract predictor variables and write to CSV.</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> CSV</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">predictors </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> collect</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">extract</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(se_aus, coords))</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">CSV</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">write</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;burramys_parvus_predictors.csv&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, predictors)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>&quot;burramys_parvus_predictors.csv&quot;</span></span></code></pre></div><p>Or convert them to a <code>DataFrame</code>:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> DataFrames</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">df </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> DataFrame</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(predictors)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">df[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">5</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">,:]</span></span></code></pre></div>`,22),e=[l];function h(k,d,r,o,c,g){return i(),a("div",null,e)}const u=s(t,[["render",h]]);export{y as __pageData,u as default};