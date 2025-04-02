import{_ as a,c as i,o as e,aA as t}from"./chunks/framework.Z4_RI1iK.js";const g=JSON.parse(`{"title":"","description":"","frontmatter":{"layout":"home","hero":{"name":"Rasters.jl","text":"Manipulating spatial data","tagline":"a powerful package that simplifies the handling of rasterized spatial data in Julia.","image":{"src":"/logo.png","alt":"Rasters"},"actions":[{"theme":"brand","text":"Get Started","link":"/get_started"},{"theme":"alt","text":"View on Github","link":"https://github.com/rafaqz/Rasters.jl"},{"theme":"alt","text":"API Reference","link":"/api"}]},"features":[{"title":"💥💠 Core Functionality","details":"Defines common types and methods for reading, writing, and manipulating rasterized spatial data. <a class=\\"highlight-link\\">Rasters.jl</a> provides unified data handling through types like <a class=\\"highlight-link\\">Raster</a>, <a class=\\"highlight-link\\">RasterStack</a>, and <a class=\\"highlight-link\\">RasterSeries</a>, offering seamless abstraction regardless of storage backend.","link":"/manual/methods"},{"title":"⚙️ Data-Source Abstraction","details":"Rasters provides a standardized interface that enables various data source types to be used with consistent syntax. The data can include <a class=\\"highlight-link\\">GeoTIFF</a> or <a class=\\"highlight-link\\">NetCDF</a> files, in-memory Arrays, or <a class=\\"highlight-link\\">CuArrays</a> on the <a class=\\"highlight-link\\">GPU</a>—all of which will behave in the same way.","link":"/manual/data_sources"},{"title":"🗂️🌐 Data Formats","details":"A <a class=\\"highlight-link\\">RasterStack</a> can be backed by a <a class=\\"highlight-link\\">NetCDF</a> or <a class=\\"highlight-link\\">HDF5</a> file, a NamedTuple of Rasters holding <a class=\\"highlight-link\\">.tif</a> files, or all Rasters in memory. Users do not need to worry about the specifics of spatial file types."},{"title":"🌍🔍 Effortless Spatial Lookups","details":"<a class=\\"highlight-link\\">Projected</a> lookups with cylindrical projections can be indexed using other cylindrical projections by setting the <a class=\\"highlight-link\\">mappedcrs</a> keyword during construction. You don’t need to know the underlying projection, as the conversion is handled automatically. This means that <a class=\\"highlight-link\\">lat/lon EPSG(4326)</a> can be used seamlessly if needed."},{"title":"🧩⚡ Modular Extensions with Faster Loading","details":"From version 0.8, Rasters.jl uses Julia's package extension system. On Julia 1.9, packages can be placed in <a class=\\"highlight-link\\">extensions</a> to load only when needed, reducing loading time to under a second. However, each backend or feature <a class=\\"highlight-link\\">must be loaded manually</a>.","link":"/#Package-extensions"},{"title":"🐞📌 Bug reports","details":"Raster data is complex, and there are many opportunities for subtle or obvious bugs to appear. We need bug reports to reduce how often these issues occur over time. However, it's also crucial that issues are easy to reproduce; otherwise, it becomes impractical to fix them. Please follow <a href=\\"#Bugs\\" class=\\"highlight-link\\">these guidelines</a>!","link":"/#Bugs"}]},"headers":[],"relativePath":"index.md","filePath":"index.md","lastUpdated":null}`),l={name:"index.md"};function n(h,s,p,o,d,r){return e(),i("div",null,s[0]||(s[0]=[t("",12)]))}const c=a(l,[["render",n]]);export{g as __pageData,c as default};
