import{T as p}from"./chunks/theme.uqQvuT1N.js";import{R as s,an as i,ao as u,ap as c,aq as l,ar as f,as as d,at as m,au as h,av as g,aw as A,d as v,u as y,v as w,s as C,ax as P,ay as R,az as T,am as b}from"./chunks/framework.Z4_RI1iK.js";function r(e){if(e.extends){const a=r(e.extends);return{...a,...e,async enhanceApp(t){a.enhanceApp&&await a.enhanceApp(t),e.enhanceApp&&await e.enhanceApp(t)}}}return e}const n=r(p),E=v({name:"VitePressApp",setup(){const{site:e,lang:a,dir:t}=y();return w(()=>{C(()=>{document.documentElement.lang=a.value,document.documentElement.dir=t.value})}),e.value.router.prefetchLinks&&P(),R(),T(),n.setup&&n.setup(),()=>b(n.Layout)}});async function S(){globalThis.__VITEPRESS__=!0;const e=D(),a=x();a.provide(u,e);const t=c(e.route);return a.provide(l,t),a.component("Content",f),a.component("ClientOnly",d),Object.defineProperties(a.config.globalProperties,{$frontmatter:{get(){return t.frontmatter.value}},$params:{get(){return t.page.value.params}}}),n.enhanceApp&&await n.enhanceApp({app:a,router:e,siteData:m}),{app:a,router:e,data:t}}function x(){return A(E)}function D(){let e=s;return h(a=>{let t=g(a),o=null;return t&&(e&&(t=t.replace(/\.js$/,".lean.js")),o=import(t)),s&&(e=!1),o},n.NotFound)}s&&S().then(({app:e,router:a,data:t})=>{a.go().then(()=>{i(a.route,t.site),e.mount("#app")})});export{S as createApp};
