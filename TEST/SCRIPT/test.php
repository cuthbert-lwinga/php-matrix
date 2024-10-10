<?php

ini_set('memory_limit', '4G');
error_reporting(E_ALL);
ini_set('display_errors', 1);

// Function to compare float values with tolerance
function floatEquals($a, $b, $tolerance = 1e-6) {
    return abs($a - $b) < $tolerance;
}

// Function to compare matrices with tolerance and return differences
function matrixDiff($matrix1, $matrix2, $tolerance = 1e-6) {
    $diff = [];
    if (!is_array($matrix1) || !is_array($matrix2)) {
        return ["Matrices are not comparable"];
    }
    if (count($matrix1) != count($matrix2) || count($matrix1[0]) != count($matrix2[0])) {
        return ["Matrix dimensions do not match"];
    }
    for ($i = 0; $i < count($matrix1); $i++) {
        for ($j = 0; $j < count($matrix1[0]); $j++) {
            if (!floatEquals($matrix1[$i][$j], $matrix2[$i][$j], $tolerance)) {
                $diff[] = "[$i][$j]: {$matrix1[$i][$j]} << {$matrix2[$i][$j]}";
            }
        }
    }
    return $diff;
}

// Function to print colored text
function coloredText($text, $color) {
    $colors = [
        'red' => "\033[0;31m",
        'green' => "\033[0;32m",
        'yellow' => "\033[0;33m",
        'blue' => "\033[0;34m",
        'orange' => "\033[0;33m",
        'reset' => "\033[0m",
    ];
    return $colors[$color] . $text . $colors['reset'];
}

// Dynamically construct the path to the JSON file
$currentDir = __DIR__;
$jsonFile = dirname($currentDir) . '/DATA/matrix_operations_results.json';

// Check if the file exists
if (!file_exists($jsonFile)) {
    die(coloredText("Error: JSON file not found at $jsonFile\n", 'red'));
}

$jsonData = json_decode(file_get_contents($jsonFile), true);

if (json_last_error() !== JSON_ERROR_NONE) {
    die(coloredText("Error: Failed to parse JSON file\n", 'red'));
}

//$jsonData = json_decode(file_get_contents($jsonFile), true);

echo coloredText("Matrix Operations Test\n", 'blue');
echo str_repeat('-', 40) . "\n";

// Initialize the Matrix object with the initial matrix
$initialMatrix = new Matrix($jsonData['initial_matrix']);



$totalOperations = 0;
$implementedOperations = 0;
$passedOperations = 0;

// Get all methods of the Matrix class
$matrixMethods = get_class_methods('Matrix');
var_dump(implode(" ",$matrixMethods));
foreach ($jsonData as $operation => $data) {
    
    if ($operation === 'initial_matrix' || $operation === 'same_shape_matrix' || $operation === 'add_matrix') {
        continue;
    }

    $totalOperations++;
    echo coloredText("Testing operation: $operation\n", 'yellow');
    
    $pythonResult = $data['result'];
    $pythonTime = $data['execution_time'];
    
    $startTime = microtime(true);
    
    $methodName = $operation;
    
    if($methodName=="multiply"){
        $methodName = "mul";
    }
    
    // Check if the method exists
    if (in_array($methodName, $matrixMethods) || in_array($operation, ['add_same_shape', 'add_row', 'subtract_same_shape', 'subtract_row'])) {

        
        try {
            switch ($operation) {
                case 'add_same_shape':
                    $phpResult = $initialMatrix->add(new Matrix($jsonData['same_shape_matrix']));
                    break;
                case 'add_row':
                    $phpResult = $initialMatrix->add(new Matrix($jsonData['add_matrix']));
                    break;
                case 'subtract_same_shape':
                    $phpResult = $initialMatrix->sub(new Matrix($jsonData['same_shape_matrix']));
                    break;
                case 'subtract_row':
                    $phpResult = $initialMatrix->sub(new Matrix($jsonData['add_matrix']));
                    break;
                case 'dot':
                    $phpResult = $initialMatrix->dot($initialMatrix);
                    break;
                case 'multiply':
                    $phpResult = $initialMatrix->mul($initialMatrix);
                    break;
                case 'div':
                    $phpResult = $initialMatrix->div($initialMatrix);
                    break;
                case 'exp':
                    $phpResult = $initialMatrix->exp();
                    break;
                case 'inverse':
                    $phpResult = $initialMatrix->inverse();
                    break;
                case 'determinant':
                    $phpResult = $initialMatrix->determinant();
                    break;
                case 'log':
                    $phpResult = $initialMatrix->log();
                    break;
                case 'eigen':
                    $phpResult = $initialMatrix->eigen();
                    $phpResult[0]->display();
                    break;
                case 'transpose':
                    $phpResult = $initialMatrix->transpose();
                    break;
                case 'argmax':
                    $phpResult = $initialMatrix->argmax();
                    break;
                case 'clip':
                    $phpResult = $initialMatrix->clip(0.3, 0.7);
                    break;
                case 'round':
                    $phpResult = $initialMatrix->round();
                    break;
                case 'zeros_like':
                    $phpResult = Matrix::zeros_like($initialMatrix);
                    break;
                case 'ones_like':
                    $phpResult = Matrix::ones_like($initialMatrix);
                    break;
                case 'relu':
                    $phpResult = $initialMatrix->relu(0,0);
                    break;
                case 'reshape':
                    $phpResult = $initialMatrix->reshape([400, 100]);
                    break;
                case 'sum':
                    $phpResult = $initialMatrix->sum();
                    break;
                case 'pow':
                    $phpResult = $initialMatrix->pow(2);
                    break;
                case 'abs':
                    $phpResult = $initialMatrix->abs();
                    break;
                case 'sqrt':
                    $phpResult = $initialMatrix->sqrt();
                    break;
                case 'eye':
                    $phpResult = Matrix::eye(200);
                    break;
                case 'mean':
                    $phpResult = $initialMatrix->mean();
                    break;
                case 'std':
                    $phpResult = $initialMatrix->std();
                    break;
                case 'sign':
                    $phpResult = $initialMatrix->sign();
                    break;
                case 'min':
                    $phpResult = $initialMatrix->min();
                    break;
                case 'max':
                    $phpResult = $initialMatrix->max();
                    break;
                case 'slice':
                    $phpResult = $initialMatrix->getSlice([100, 150]);
                    break;
                case 'rank':
                    $phpResult = $initialMatrix->rank();
                    break;
                case 'trace':
                    $phpResult = $initialMatrix->trace();
                    break;
                case 'frobenius_norm':
                    $phpResult = $initialMatrix->norm('fro');
                    break;
                case 'l2_norm':
                    $phpResult = $initialMatrix->norm();
                    break;
                case 'infinity_norm':
                    $phpResult = $initialMatrix->norm(INF);
                    break;
                case 'hadamard_product':
                    $phpResult = $initialMatrix->mul(new Matrix($jsonData['same_shape_matrix']));
                    break;
                case 'kronecker_product':
                    $phpResult = $initialMatrix->kron(new Matrix($jsonData['same_shape_matrix']));
                    break;
                case 'matrix_power':
                    $phpResult = $initialMatrix->pow(3);
                    break;
                case 'lu_decomposition':
                    $phpResult = $initialMatrix->lu();
                    break;
                case 'qr_decomposition':
                    $phpResult = $initialMatrix->qr();
                    break;
                case 'cholesky_decomposition':
                    $phpResult = $initialMatrix->cholesky();
                    break;
                case 'svd':
                    $phpResult = $initialMatrix->svd();
                    break;
                case 'pseudo_inverse':
                    $phpResult = $initialMatrix->pinv();
                    break;
                case 'matrix_exp':
                    $phpResult = $initialMatrix->expm();
                    break;
                case 'matrix_log':
                    $phpResult = $initialMatrix->logm();
                    break;
                case 'diagonalization':
                    $phpResult = $initialMatrix->eig();
                    break;
                case 'flatten':
                    $phpResult = $initialMatrix->flatten();
                    break;
                case 'pad':
                    $phpResult = $initialMatrix->pad(1);
                    break;
                case 'random_permute':
                    $phpResult = $initialMatrix->permute();
                    break;
                case 'sort':
                    $phpResult = $initialMatrix->sort();
                    break;
                case 'symmetry_check':
                    $phpResult = $initialMatrix->isSymmetric();
                    break;
                case 'covariance':
                    $phpResult = $initialMatrix->cov();
                    break;
                case 'correlation':
                    $phpResult = $initialMatrix->corr();
                    break;
                case 'moving_average':
                    $phpResult = $initialMatrix->movingAverage(5);
                    break;
                default:
                $phpResult = null;
                echo coloredText("Result: NOT IMPLEMENTED\n", 'orange');
                    break;
            }
            
            $implementedOperations++;
            $endTime = microtime(true);
            $phpTime = $endTime - $startTime;
            
            // Compare results
            $phpResultData = is_object($phpResult) && method_exists($phpResult, 'getData') ? $phpResult->getData() : $phpResult;
            $diff = is_array($pythonResult) ? matrixDiff($phpResultData, $pythonResult) : (!floatEquals($phpResultData, $pythonResult) ? ["Value: $phpResultData << $pythonResult"] : []);
            
            if (empty($diff)) {
                echo coloredText("Result: PASS\n", 'green');
                $passedOperations++;
            } else {
                echo coloredText("Result: FAIL\n", 'red');
                echo "Differences:\n";
                foreach (array_slice($diff, 0, 5) as $d) {
                    echo "  $d\n";
                }
                if (count($diff) > 5) {
                    echo "  ... and " . (count($diff) - 5) . " more differences\n";
                }
            }
            
            echo "Python time: " . number_format($pythonTime, 6) . " seconds\n";
            echo "PHP time: " . number_format($phpTime, 6) . " seconds\n";
            echo "Time difference: " . number_format(abs($phpTime - $pythonTime), 6) . " seconds\n";
            
            if ($phpTime < $pythonTime) {
                echo coloredText("PHP is faster\n", 'green');
            } else {
                echo coloredText("PHP is slower\n", 'red');
            }


        } catch (Throwable $e) {
            echo coloredText("Result: ERROR\n", 'red');
            echo "Error: " . $e->getMessage() . "\n";
        }
    } else {
        echo coloredText("Result: NOT IMPLEMENTED $operation \n", 'orange');
    }
    
    echo str_repeat('-', 40) . "\n";
}

echo coloredText("Matrix Operations Test Completed\n", 'blue');
echo coloredText("Summary:\n", 'yellow');
echo "Total operations: $totalOperations\n";
echo "Implemented operations: $implementedOperations\n";
echo coloredText("Passed operations: $passedOperations\n", 'green');
echo coloredText("Failed operations: " . ($implementedOperations - $passedOperations) . "\n", 'red');
echo "Pass rate: " . ($implementedOperations > 0 ? number_format(($passedOperations / $implementedOperations) * 100, 2) : 0) . "%\n";

?>