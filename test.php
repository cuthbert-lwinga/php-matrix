<?php
// Start the timer
// Generate random matrices
$rows = 5000;
$cols = 5000;


$m1 = array_fill(0, $rows, array_fill(0, $cols, mt_rand() / mt_getrandmax()));
$m2 = array_fill(0, $rows, array_fill(0, $cols, mt_rand() / mt_getrandmax()));

$start_time = microtime(true);

$matrix1 = new MatrixWrapper($m1);
$matrix2 = new MatrixWrapper($m2);

$end_time = microtime(true);
$setup_time = $end_time - $start_time;

// Output execution times
echo "Setup Time: " . $setup_time . " seconds\n";


$start_time = microtime(true);
$resultSubtract = $matrix1->sub($matrix2);
$end_time = microtime(true);
$subtract_time = $end_time - $start_time;

echo "Subtraction Time: " . $subtract_time . " seconds\n";

// Perform operations
$resultAdd = $matrix1->add($matrix2);
$end_time = microtime(true);
$add_time = $end_time - $start_time;

echo "Addition Time: " . $add_time . " seconds\n";

$resultAdd = $matrix1->add(2);
$end_time = microtime(true);
$add_time = $end_time - $start_time;

echo "Addition constant Time: " . $add_time . " seconds\n";


$resultAdd = $matrix1->sub(2);
$end_time = microtime(true);
$add_time = $end_time - $start_time;

echo "Subtraction constant Time: " . $add_time . " seconds\n";


$start_time = microtime(true);
$resultDot = $matrix1->dot($matrix2);
$end_time = microtime(true);
$dot_time = $end_time - $start_time;

echo "Dot Product Time: " . $dot_time . " seconds\n";
?>
