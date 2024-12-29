import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
}
const navTemp = {
  nav: [
    { text: 'Home', link: '/' },
    { text: 'Get Started', link: '/get_started' },
    { text: 'Manual',
      items: [
        { text: 'Methods', link: '/manual/methods' },
        { text: 'cellarea', link: '/manual/cellarea' },
        { text: 'Array Operations', link: '/manual/array_operations' },
        { text: 'Data Sources', link: '/manual/data_sources' },
        { text: 'Plots',
          items: [
            { text: 'Plots.jl', link: '/manual/plotting' },
            { text: 'Makie.jl', link: 'manual/plot_makie' },
          ]
        },  
      ]
    },
    { text: 'Tutorials',
      items: [
        { text: 'Spatial mean', link: '/tutorials/spatial_mean' },
        { text: 'Reprojection and resampling', link: '/tutorials/resample_warp'},
        { text: 'Species Distribution Modelling', link: '/tutorials/gbif_wflow' },
      ]
    },
    { text: 'Ecosystem',
      items: [
        { text: 'DimensionalData.jl', link: 'https://rafaqz.github.io/DimensionalData.jl' },
        { text: 'DiskArrays.jl', link: 'https://github.com/JuliaIO/DiskArrays.jl' },
        { text: 'GeoInterface.jl', link: 'https://github.com/JuliaGeo/GeoInterface.jl' },
        { text: 'NCDatasets.jl', link: 'https://alexander-barth.github.io/NCDatasets.jl/stable/' },
        { text: 'ArchGDAL.jl', link: 'https://github.com/yeesian/ArchGDAL.jl' },
        { text: 'GRIBDatasets.jl', link: 'https://github.com/JuliaGeo/GRIBDatasets.jl' },
        { text: 'ZarrDatasets.jl', link: 'https://github.com/JuliaGeo/ZarrDatasets.jl' },
      ]
    },
    { text: 'API', link: '/api' }
  ],
}

const nav = [
  ...navTemp.nav,
  {
    component: 'VersionPicker'
  }
]

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  title: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  description: "Manipulating rasterized spatial data",
  lastUpdated: true,
  cleanUrls: true,
  outDir: 'REPLACE_ME_DOCUMENTER_VITEPRESS', // This is required for MarkdownVitepress to work correctly...
  head: [
    ['link', { rel: 'icon', href: '/favicon.ico' }],
    ['script', {src: `${getBaseRepository(baseTemp.base)}versions.js`}],
    ['script', {src: `${baseTemp.base}siteinfo.js`}]
  ],
  ignoreDeadLinks: true,

  markdown: {
    math: true,
    config(md) {
      md.use(tabsMarkdownPlugin),
      md.use(mathjax3),
      md.use(footnote)
    },
    theme: {
      light: "github-light",
      dark: "github-dark"}
  },

  themeConfig: {
    outline: 'deep',
    // https://vitepress.dev/reference/default-theme-config
    logo: '/logo.png',
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar: [
      { text: 'Get Started', link: '/get_started' },
      { text: 'Manual',
        items: [
          { text: 'Methods', link: '/manual/methods' },
          { text: 'cellarea', link: '/manual/cellarea' },
          { text: 'Array Operations', link: '/manual/array_operations' },
          { text: 'Data Sources', link: '/manual/data_sources' },
          { text: 'Plots',
            items: [
              { text: 'Plots.jl', link: '/manual/plotting' },
              { text: 'Makie.jl', link: '/manual/plot_makie' },
            ]
          },  
        ]
      },
      { text: 'Tutorials',
        items: [
          { text: 'Spatial mean', link: '/tutorials/spatial_mean' },
          { text: 'Reprojection and resampling', link: '/tutorials/resample_warp'},
          { text: 'Species Distribution Modelling', link: '/tutorials/gbif_wflow' },
        ]
      },
      { text: 'API', link: '/api' }
    ],
    editLink: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    socialLinks: [
      // { icon: 'github', link: 'REPLACE_ME_DOCUMENTER_VITEPRESS' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/" target="_blank"><strong>DocumenterVitepress.jl</strong></a> <br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}. Released under the MIT License.`
    }
  }
})
