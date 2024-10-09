import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  title: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  description: "Manipulating rasterized spatial data",
  lastUpdated: true,
  cleanUrls: true,
  outDir: 'REPLACE_ME_DOCUMENTER_VITEPRESS', // This is required for MarkdownVitepress to work correctly...
  head: [['link', { rel: 'icon', href: '/public/favicon.ico' }]],
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
    nav: [
      { text: 'Home', link: '/' },
      { text: 'Get Started', link: '/get_started' },
      { text: 'Methods',
        items: [
          { text: 'Overview', link: '/methods' },
          { text: 'Array Operations', link: '/array_operations' },
        ]
       },
       { text: 'Data Sources',
       items: [
         { text: 'Overview', link: '/data_sources' },
         { text: 'GBIF', link: '/gbif_wflow' }
        ]
      },
      { text: 'Plots',
      items: [
        { text: 'Plots.jl', link: '/plotting' },
        { text: 'Makie.jl', link: '/plot_makie' },
      ]
     },
     { text: 'Ecosystem',
      items: [
        { text: 'DimensionalData.jl', link: 'https://rafaqz.github.io/DimensionalData.jl/dev/' },
        { text: 'NCDatasets.jl', link: 'https://alexander-barth.github.io/NCDatasets.jl/stable/' },
        { text: 'ArchGDAL.jl', link: 'https://yeesian.com/ArchGDAL.jl/stable/' },
        { text: 'HDF5.jl', link: 'https://juliaio.github.io/HDF5.jl/stable/' },
       ]
     },
      { text: 'API', link: '/api' }
    ],

    sidebar: [
      { text: 'Get Started', link: '/get_started' },
      { text: 'Methods',
        items: [
          { text: 'Overview', link: '/methods' },
          { text: 'Array Operations', link: '/array_operations' },
        ]
       },
       { text: 'Data Sources',
       items: [
         { text: 'Overview', link: '/data_sources' },
         { text: 'GBIF', link: '/gbif_wflow' }
        ]
      },
      { text: 'Plots',
      items: [
        { text: 'Plots.jl', link: '/plotting' },
        { text: 'Makie.jl', link: '/plot_makie' },
      ]
     },

      { text: 'API', link: '/api' }
    ],
    editLink: {
      pattern: 'REPLACE_ME_DOCUMENTER_VITEPRESS'
    },
    socialLinks: [
      { icon: 'github', link: 'REPLACE_ME_DOCUMENTER_VITEPRESS' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/" target="_blank"><strong>DocumenterVitepress.jl</strong></a> <br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}. Released under the MIT License.`
    }
  }
})