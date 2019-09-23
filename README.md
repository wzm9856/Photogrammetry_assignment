# Photogrammetry_project

是摄影测量后方交会实验用程序，总共包括三个文件：

[source.cpp](https://github.com/wzm9856/Photogrammetry_project/blob/master/source.cpp)为程序主文件，包括后方交会的主要算法。题目用的四个点位数据在main函数中作为四个自定的Point类的对象被声明，用此方法可方便调用和更改数据。

[matrix.h](https://github.com/wzm9856/Photogrammetry_project/blob/master/matrix.h)和[matrix.cpp](https://github.com/wzm9856/Photogrammetry_project/blob/master/matrix.cpp)为主文件使用的矩阵类。由于c++中并无自带的矩阵运算库，所以自建了这个Matrix类，其中包括矩阵的加减法，乘法，打印矩阵和矩阵求逆等方法，全部为自己手打且经过自己验证可用的。
