# PHP Matrix Operations Library 🧮🚀

This PHP library provides a high-performance implementation of matrix operations, designed for efficient handling of large matrices. It utilizes parallel processing to speed up computations and includes a wide range of matrix operations.

## Features 🌟

- Creation of random matrices 🎲
- Basic arithmetic operations (addition, subtraction, multiplication, division) ➕➖✖️➗
- Advanced matrix operations (transpose, inverse, determinant, eigenvalues) 🔄🔢
- Element-wise operations (exp, log, clip, relu) 📊
- Statistical functions (mean, sum, standard deviation) 📈
- Shape manipulation (reshape, slicing) ✂️🔁
- Initialization functions (zeros, ones, eye, glorot_uniform) 🔢
- Performance optimizations for large matrices 🏎️💨

## Requirements 📋

- PHP 7.4 or higher 🐘
- PCNTL extension enabled 🔌
- Large amount of memory (adjustable via `memory_limit` setting) 💾
- PHP-CPP library 🔧

## Usage 🛠️

This package supports operation chaining, making it possible to do operations like `->pow()->dot()` ... and so forth 🔗

Here's a basic example of how to use the library:

```php
<?php
// Create a random matrix
$matrix = new Matrix();
$matrix->random(1000, 1000, 0.2, 0.8);

// Perform operations
$result = $matrix->transpose()->dot($matrix);

// Display results
$result->display();
```

## PHP-CPP Usage 🔧🐘

This library leverages PHP-CPP to achieve high performance. Here's how to use PHP-CPP with this library:

1. Install PHP-CPP:
   ```
   git clone https://github.com/CopernicaMarketingSoftware/PHP-CPP.git
   cd PHP-CPP
   make
   sudo make install
   ```

2. Compile the extension:
   ```
   g++ -std=c++11 -Wall -c -O2 -fpic -o matrix.o matrix.cpp
   g++ -shared -o matrix.so matrix.o -lphpcpp
   ```

3. Add the extension to your php.ini file:
   ```
   extension=matrix.so
   ```

4. Restart your PHP server to load the new extension.

Now you can use the Matrix class in your PHP code as shown in the usage example above! 🎉

## Performance Considerations 🏎️

- The library uses parallel processing to handle large matrices efficiently. 🚀
- Memory usage can be high for large matrices. Adjust the `memory_limit` setting as needed. 💾
- Performance metrics are provided for various operations to help optimize your code. ⏱️

## Available Operations 🧰

- Matrix creation: `random()`, `zeros()`, `ones()`, `eye()`, `glorot_uniform()` 🆕
- Arithmetic: `add()`, `sub()`, `mul()`, `div()`, `dot()` 🧮
- Element-wise: `exp()`, `log()`, `clip()`, `relu()`, `abs()`, `sqrt()`, `pow()` 📊
- Statistical: `mean()`, `sum()`, `std()`, `min()`, `max()` 📈
- Shape: `reshape()`, `transpose()`, `getSlice()`, `getRow()`, `getCol()`, `setRow()`, `setCol()` 🔄
- Linear algebra: `inverse()`, `determinant()`, `eigen()` 🔢
- Utility: `shape()`, `display()`, `getData()` 🛠️

## Contributing 🤝

Contributions to improve the library are welcome. Please ensure that any pull requests include appropriate tests and documentation. Let's make this library awesome together! 💪

## License 📄

Happy matrix computing! 🎉🧮🚀