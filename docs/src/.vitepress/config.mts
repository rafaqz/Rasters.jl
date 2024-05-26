import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: '/Rasters.jl/',
  title: "Rasters",
  description: "Manipulating rasterized spatial data",
  lastUpdated: true,
  cleanUrls: true,
  outDir: '../final_site', // This is required for MarkdownVitepress to work correctly...
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
    logo: { src: '/logo.png', width: 24, height: 24 },
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
      pattern: 'https://github.com/MakieOrg/Tyler.jl/edit/master/docs/src/:path'
    },
    socialLinks: [
      { icon: 'github', link: 'https://github.com/MakieOrg/Tyler.jl' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/" target="_blank"><strong>DocumenterVitepress.jl</strong></a> <br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}. Released under the MIT License.`
    }
  }
})