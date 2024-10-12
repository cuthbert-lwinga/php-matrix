<?php
ini_set('memory_limit', '400G'); // Increase to 2GB or more as needed
error_reporting(E_ALL);
ini_set('display_errors', 1);

// Function to measure execution time of a given operation and return the time
function measure_time($callback) {
    $start_time = microtime(true);
    $callback();
    $end_time = microtime(true);
    return $end_time - $start_time;
}

function test_dot_product_performance($rows, $cols, $min_val, $max_val, $outputFile) {
    // Initialize matrices using random method
    $matrix1 = new Matrix();
    $matrix1->random($rows, $cols, $min_val, $max_val);
    
    $matrix2 = new Matrix();
    $matrix2->random($rows, $cols, $min_val, $max_val);

    // Measure dot product time
    $dot_product_time = measure_time(function() use ($matrix1, $matrix2) {
        $matrix1->dot($matrix2);
    });

    // Save time is set to 0 as there's no specific save operation
    $save_time = 0;

    // Check if the file exists to write the header
    $file_exists = file_exists($outputFile);
    $file = fopen($outputFile, 'a');
    if ($file !== false) {
        if (!$file_exists) {
            // Write the CSV header if the file does not exist
            fputcsv($file, ['Matrix Size', 'Save Time', 'Dot Product Time']);
        }
        echo "\t$rows,\t$save_time,\t$dot_product_time \n";
        // Append results to the CSV file
        fputcsv($file, [$rows, $save_time, $dot_product_time]);
        fclose($file);
    } else {
        echo "Failed to open the output file for writing.\n";
    }
}


$min_val = 0.2;
$max_val = 0.8;

for($j = 1; $j <= 10;$j++){
    $scaling = $j*0.1;

    echo "\n SCALLING($scaling) \n";

    $outputFile = "performance/matrix_performance_with_php_cpp_scalling_$scaling.csv";
    Matrix::setThreadScalingFactor($scaling);
    echo "'Matrix Size',\t'Save Time',\t'Dot Product Time'\n";

    for($i = 1; $i < 20000;$i+=10){

        test_dot_product_performance($i, $i, $min_val, $max_val, $outputFile);

    }
}
?>
