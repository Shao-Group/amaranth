This folder contains a small example data `example-input.bam` for users to test with Amaranth. The example data is also available in [Release](https://github.com/Shao-Group/amaranth/releases/tag/v0.1.0).

Please use the following link to download from GitHub Release:

```bash
wget https://github.com/Shao-Group/amaranth/releases/download/v0.1.0/example-input.bam
```

> ‼️ Note that `wget` works only with GitHub Release, but it does NOT work with GitHub repo files. It will be truncated if you use wget to download GitHub files directly from a repo,

You can use `sha256` to test the integrity of downloaded data.

```bash
sha256 example-input.bam
```

The output should be `59e036720e6539336600409bb7a466dd82a51e8ad98a30016951aea42f21fba6`