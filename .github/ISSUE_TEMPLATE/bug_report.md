---
name: Bug report
about: Describe your issue here, including a minimum working example all files.
title: ''
labels: bug
assignees: ''

---

Raster data is complicated, there are a lot of formats and a lot of things that can go wrong. So we can't usually fix your problem unless you help us by filling out the code blocks below. 

Its important that your example just runs, and we can copy and paste the code into Julia and it does that same thing for us that it does for you. That means either generating or downloading data inside the script, e.g. with `download`. 

The faster we can get to your problem, the more likely it will be fixed quickly. If we have 15 minutes of free dev time for your issue, you want all 15 to be looking at code, not 10 trying to download some file from a website and making sure the file path matches in the script, leaving 5 minutes to look at the error and fix it.

```julia
All the julia code needed to make the error goes here including download or generation of all data used

This is the most important thing in the issue, we cant fix bugs without it.
```

```julia
The full complete error from start to finish goes here

This is the second most important thing in the issue, we can't be sure we fixed the same problem unless we can see this.
```

```julia
Paste the results of `versioninfo()` here.

Paste the results of `] status Rasters` here.
```
