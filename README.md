# usp654-f17

Source code for course website for USP654 at Portland State University, Fall 2017.

Course website is at https://usp654.github.io

# Steps seting up a new course website

1. Fork the repo at https://github.com/usp654/usp654_website;
1. Rename the forked repo to whatever you need, for example, https://github.com/<course#>/<course#>_website;
1. Create a new repo on github <course#>.github.io (note that you have to own https://github.com/<course#> for this to work);
1. Clone https://github.com/<course#>/<course#>_website.git as a local repo; rename usp654_website.Rproj to <course#>_website.Rproj;
1. Add <course#>.github.io as a submodule: `git submodule add git@github.com:<course#>/<course#>.github.io.git docs`;
1. Customize the website by modifying `config.toml`, files in `content` and `static`;
1. Build, push changes and publish the website with `make build`, `make push`, and `make publish` respectively. If you want better commit messages customize the corresponding commands in `Makefile`

# Attributions

- This website is adapted from [STA 112FS by Mine Çetinkaya-Rundel](http://www2.stat.duke.edu/courses/Fall17/sta112.01/), which in turn is based on [ESPM-157 by Carl Boettinger](https://espm-157.carlboettiger.info/).
- Background photo [~Portland at Dusk~ by Alejandro Rdguez](https://www.flickr.com/photos/90642235@N04/19400834865/)# stat564.github.io
# stat564_website
# stat564_website
