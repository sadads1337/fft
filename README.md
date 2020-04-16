FFT Bachelor degree work
==============

Something...

Dependencies
-----

* **C++17**.
* [**Python 2.7**](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=2ahUKEwiTtuuZrNboAhVcwsQBHWwcClkQFjAAegQIAhAB&url=https%3A%2F%2Fwww.python.org%2Fdownload%2Freleases%2F2.7%2F&usg=AOvVaw0zlkFAAPj_mtiLEq6iCUxh) + matplotlib + numpy. Download then call for `pip2 install numpy matplotlib`
* [**Intel MKL**](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=2ahUKEwiwv-ncrNboAhXkx6YKHSDUDVgQFjAAegQIFRAC&url=https%3A%2F%2Fsoftware.intel.com%2Fen-us%2Fmkl&usg=AOvVaw2E4mKkrupU-h-MctqNhKyu) must be installed manually.
* [**GTest**](https://github.com/google/googletest) will be downloaded automatically.
* [**matplotlib-cpp**](https://github.com/lava/matplotlib-cpp) forked as submodule.
* [**benchmark**](https://github.com/google/benchmark) as submodule.
* **OpenMP** optional.

Build
-----
**Minimal** build example via **console** on MacOS + apple-clang.
```c++ based
git clone <url-for-this-repo>
git submodule update
mkdir build
cd build
cmake .. -G Ninja
ninja -j 16
```
