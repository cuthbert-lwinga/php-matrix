<?php
ini_set('memory_limit', '200G'); // Increase to 2GB or more as needed

// Helper function to create a random matrix
function create_random_matrix($rows, $cols) {
    $matrix = array_fill(0, $rows, array_fill(0, $cols, 0));
    for ($i = 0; $i < $rows; ++$i) {
        for ($j = 0; $j < $cols; ++$j) {
            $matrix[$i][$j] = mt_rand() / mt_getrandmax();
        }
    }
    return $matrix;
}

// Function to measure execution time of a given operation
function measure_time($callback, $description) {
    $start_time = microtime(true);
    $callback();
    $end_time = microtime(true);
    $execution_time = $end_time - $start_time;
    echo $description . ": " . $execution_time . " seconds\n";
}

// Generate random matrices
$rows = 10000; // Adjust size as needed
$cols = 10000;
$min_val = 0.2;
$max_val = 0.8;

echo " MATRIX($rows,$cols)\n\n";

$m1 = create_random_matrix($rows, $cols);
$m2 = create_random_matrix($rows, $cols);

$matrix1 = new MatrixWrapper($m1);
$matrix2 = new MatrixWrapper($m2);

measure_time(function() use ($m1) {
    $matrix1 = new MatrixWrapper($m1);
}, "Init Time");


measure_time(function() use ($matrix1) {
    $matrix1->shape();
}, "Finding Shape Time");

// Perform subtraction first and then addition
measure_time(function() use ($matrix1, $matrix2) {
    $matrix1->sub($matrix2);
}, "Subtraction Time");

measure_time(function() use ($matrix1, $matrix2) {
    $matrix1->add($matrix2);
}, "Addition Time after Subtraction");

// Perform addition first and then subtraction
measure_time(function() use ($matrix1, $matrix2) {
    $matrix1->add($matrix2);
}, "Addition Time");

measure_time(function() use ($matrix1, $matrix2) {
    $matrix1->sub($matrix2);
}, "Subtraction Time after Addition");

// Measure constant addition time
measure_time(function() use ($matrix1) {
    $matrix1->add(2);
}, "Addition Constant Time");

// Measure dot product time
measure_time(function() use ($matrix1, $matrix2) {
    $matrix1->dot($matrix2);
}, "Dot Product Time");

measure_time(function() use ($matrix1, $matrix2) {
    $matrix1->div($matrix2);
}, "Matrix-to-Matrix Division Time");

// Measure matrix-to-scalar division time
measure_time(function() use ($matrix1) {
    $matrix1->div(2);
}, "Matrix-to-Scalar Division Time");

// Measure transpose time
measure_time(function() use ($matrix1) {
    $matrix1->transpose();
}, "Transpose Time");

// Measure Multiplication time
measure_time(function() use ($matrix1) {
    $matrix1->mul(2);
}, "Multiplication Time");

// Measure argmax time
measure_time(function() use ($matrix1) {
    $matrix1->argmax(1);
}, "Arg max Time");

// Measure clip time
measure_time(function() use ($matrix1, $min_val, $max_val) {
    $matrix1->clip($min_val, $max_val);
}, "Clip Time");

// Measure log
measure_time(function() use ($matrix1) {
    $matrix1->log();
}, "Log Time");

// Measure exp
measure_time(function() use ($matrix1) {
    $matrix1->exp();
}, "Exp Time");

// Measure sum
measure_time(function() use ($matrix1) {
    $matrix1->sum();
}, "Sum Time");

measure_time(function() use ($matrix1) {
    $matrix1->sum(0);
}, "Sum along columns Time");

measure_time(function() use ($matrix1) {
    $matrix1->sum(1);
}, "Sum along rows Time");

// Measure inverse
measure_time(function() use ($matrix1) {
    $matrix1->inverse();
}, "Inverse Time");

// Measure determinant
measure_time(function() use ($matrix1) {
    $matrix1->determinant();
}, "Determinant Time");

// Measure eigen
measure_time(function() use ($matrix1) {
    $matrix1->eigen();
}, "Eigenvalues and Eigenvectors Time");
?>