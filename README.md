# Statistic

Source code for course website for STAT651-6 at Portland State University, Fall 2018.

Course website is at https://stat564.netlify.com/

# Steps seting up a new course website

1. Fork the repo at https://github.com/usp654/usp654_website;
1. Rename the forked repo to whatever you need, for example, https://github.com/<course#>/<course#>_website;
1. Create a new repo on github <course#>.github.io (note that you have to own https://github.com/<course#> for this to work);
1. Clone https://github.com/<course#>/<course#>_website.git as a local repo; rename usp654_website.Rproj to <course#>_website.Rproj;
1. Add <course#>.github.io as a submodule: `git submodule add git@github.com:<course#>/<course#>.github.io.git docs`;
1. Customize the website by modifying `config.toml`, files in `content` and `static`;
1. Build, push changes and publish the website with `make build`, `make push`, and `make publish` respectively. If you want better commit messages customize the corresponding commands in `Makefile`


