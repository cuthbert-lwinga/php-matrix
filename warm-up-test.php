<?php
ini_set('memory_limit', '400G'); // Increase to 2GB or more as needed
error_reporting(E_ALL);
ini_set('display_errors', 1);
// include_once("PHPThreading/ThreadManager.php");
// Helper function to create a random matrix
// function create_random_matrix($rows, $cols) {
//     $matrix = array_fill(0, $rows * $cols, 0.0);
//     $matrixPtr = &$matrix;
//     $numElements = $rows * $cols;
    
//     for ($i = 0; $i < $numElements; $i++) {
//         $matrixPtr[$i] = mt_rand() / mt_getrandmax();
//     }
    
//     $result = array_chunk($matrix, $cols);
    
//     return $result;
// }
function create_random_matrix($rows, $cols) {
    $numProcesses = 100; // Adjust this value based on your system's available cores
    $colsPerProcess = ceil($cols / $numProcesses);
    $tempFile = tempnam(sys_get_temp_dir(), 'matrix_');
    echo "\n saving to $tempFile\n";

    for ($i = 0; $i < $numProcesses; $i++) {
        $pid = pcntl_fork();
        if ($pid == -1) {
            die("Failed to fork");
        } else if ($pid == 0) {
            // Child process
            $startCol = $i * $colsPerProcess;
            $endCol = min(($i + 1) * $colsPerProcess, $cols);
            generate_columns($startCol, $endCol, $rows, $tempFile);
            exit(0);
        }
    }

    // Parent process
    while (pcntl_waitpid(0, $status) != -1) {
        // Wait for all child processes to finish
    }

    $matrix = [];
    $file = fopen($tempFile, 'r');
    if ($file !== false) {
        while (($line = fgets($file)) !== false) {
            $matrix[] = explode(',', rtrim($line));
        }
        fclose($file);
    }

    unlink($tempFile);
    return $matrix;
}

function generate_columns($startCol, $endCol, $rows, $tempFile) {
    $file = fopen($tempFile, 'a');
    if ($file === false) {
        throw new Exception("Failed to open temporary file");
    }

    for ($col = $startCol; $col < $endCol; $col++) {
        $row = '';
        for ($row_idx = 0; $row_idx < $rows; $row_idx++) {
            $row .= (($row_idx > 0) ? ',' : '') . mt_rand() / mt_getrandmax();
        }
        fwrite($file, $row . PHP_EOL);
    }

    fclose($file);
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
$rows = 100; // Adjust size as needed
$cols = 100;
$min_val = 0.2;
$max_val = 0.8;

$m1 = create_random_matrix($rows, $cols);
//$m2 = create_random_matrix($rows, $cols);



$rows = 1000; // Adjust size as needed
$cols = 1000;
$min_val = 0.2;
$max_val = 0.8;
$scalling = 1.0;

MatrixWrapper::setThreadScalingFactor($scalling);

echo " MATRIX($rows,$cols) scalling: $scalling\n\n";

$matrix1 = new MatrixWrapper();
$matrix1->random($rows,$cols,$min_val,$max_val);

$matrix2 = $matrix1;//new MatrixWrapper($m2);

// measure_time(function() use ($m1) {
//     $matrix1 = new MatrixWrapper($m1);
// }, "Init Time");

measure_time(function() use ($m1,$rows,$cols,$min_val,$max_val) {
    $matrix1 = new MatrixWrapper();
    $matrix1->random($rows,$cols,$min_val,$max_val);
}, "Init Time for rand gen ($rows,$cols,$min_val,$max_val)");


measure_time(function() use ($matrix1) {
    ($matrix1->shape());
}, "Finding Shape Time");

// Perform subtraction first and then addition
// measure_time(function() use ($matrix1, $matrix2) {
//     $matrix1->sub($matrix2);
// }, "Subtraction Time");

// measure_time(function() use ($matrix1, $matrix2) {
//     $matrix1->add($matrix2);
// }, "Addition Time after Subtraction");

// // Perform addition first and then subtraction
// measure_time(function() use ($matrix1, $matrix2) {
//     $matrix1->add($matrix2);
// }, "Addition Time");

// measure_time(function() use ($matrix1, $matrix2) {
//     $matrix1->sub($matrix2);
// }, "Subtraction Time after Addition");

// // Measure constant addition time
// measure_time(function() use ($matrix1) {
//     $matrix1->add(2);
// }, "Addition Constant Time");

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