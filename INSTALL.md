<!--
# Install from bioconda
The easiest way to install `amaranth` is using bioconda:
```bash
conda install bioconda::amaranth-assembler`
```



# Download pre-compiled binary

We provide pre-compiled binary executable of `amaranth` (for Linux only). Please check out GitHub Releases (`amaranth-{version}-linux` in [Latest Release](https://github.com/Shao-Group/amaranth/releases/latest)).

Use `chmod` to make it executable, if not yet. You can test whether it works by printing the help message.

```bash
# make executable
chmod +x ./amaranth-*-linux
# test 
./amaranth-*-linux --help
```
-->


# Install from source 

Download the source code (`amaranth-{version}.tar.gz`) of Amaranth from [Releases](https://github.com/Shao-Group/amaranth/releases).

```bash
# For example, to download v0.1.0 source code
wget https://github.com/Shao-Group/amaranth/releases/latest/download/amaranth-0.1.0.tar.gz
tar -zxf amaranth-0.1.0.tar.gz
cd amaranth-0.1.0
```

Amaranth is implemented in C++. A compiler that supports C++20 is needed.

Amaranth uses additional libraries of Boost and htslib. 
If they have not been installed in your system, you first
need to download and install them. You might also need to
export the runtime library path to certain environmental
variable (for example, `LD_LIBRARY_PATH`, for most linux distributions).
After install these dependencies, you then compile the source code of Amaranth.
If some of the above dependencies are not installed to the default system 
directories (for example, `/usr/local`, for most linux distributions),
their corresponding installing paths should be specified to `configure` of Amaranth.
The entire installation typically takes a few minutes to complete.

## Download Boost

If Boost has not been downloaded/installed, download Boost
[(license)](http://www.boost.org/LICENSE_1_0.txt) from (http://www.boost.org).
Uncompress it somewhere (compiling and installing are not necessary).

## Install htslib

If htslib has not been installed, download htslib 
[(license)](https://github.com/samtools/htslib/blob/develop/LICENSE)
from (http://www.htslib.org/) with version 1.5 or higher.
Note that htslib relies on zlib. So if zlib has not been installed in your system,
you need to install zlib first. To do so, download zlib
[(license)](https://zlib.net/zlib_license.html) at (https://zlib.net/).
Use the following commands to install zlib:

```
./configure
make
make install
```

After installing zlib, use the following commands to build htslib:

```
./configure --disable-bz2 --disable-lzma --disable-gcs --disable-s3 --enable-libcurl=no
make
make install
```

The default installation location of htslib is `/usr/lib`.
If you would install it to a different location, replace the above `configure` line with
the following (by adding `--prefix=/path/to/your/htslib` to the end):

```
./configure --disable-bz2 --disable-lzma --disable-gcs --disable-s3 --enable-libcurl=no --prefix=/path/to/your/htslib
```

In this case, you also need to export the runtime library path (note that there
is an additional `lib` following the installation path):

```
export LD_LIBRARY_PATH=/path/to/your/htslib/lib:$LD_LIBRARY_PATH
```

## Build Amaranth

Use the following to compile:

```
./configure --with-htslib=/path/to/your/htslib --with-boost=/path/to/your/boost
make
```

If some of the dependencies are installed in the default system directory (for example, `/usr/lib`),
then the corresponding `--with-` option might not be necessary.
The executable file `amaranth` will appear at `src/amaranth`.

You can test whether amaranth is correctly installed by printing the help message:
```bash
./src/amaranth --help
```
