<?php
ini_set('memory_limit', '400G'); // Increase to 2GB or more as needed
error_reporting(E_ALL);
ini_set('display_errors', 1);



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
$rows = 1000; // Adjust size as needed
$cols = 1000;
$min_val = 0.2;
$max_val = 0.8;
$scalling = 1;
$m1 = create_random_matrix(20, 20);
$m2 = [];//$m1;//create_random_matrix($rows, $cols);

Matrix::setThreadScalingFactor($scalling);

echo " MATRIX($rows,$cols) scalling: $scalling\n\n";

$matrix1 = new Matrix($m1);
$matrix1->random($rows,$cols,$min_val,$max_val);
//$matrix1 = $matrix1->random($rows,$cols,$min_val,$max_val,$binomial = true);

$matrix2 = $matrix1;//new MatrixWrapper($m2);

measure_time(function() use ($m1) {
    echo "NOT TESTING METHODE IS SLOW FOR BIG MATRIX \n";
    //$matrix1 = new Matrix($m1);
}, "Init Time with passed array: This should be avoided for massive matrix it's computationally intensive");

measure_time(function() use ($m1,$rows,$cols,$min_val,$max_val) {
    $matrix1 = new Matrix();
    $matrix1->random($rows,$cols,$min_val,$max_val);
}, "Init Time for rand gen ($rows,$cols,$min_val,$max_val)");


measure_time(function(){
    var_dump((new Matrix([1,2,3,4,5,6]))->shape());
}, "shape from Matrix");

measure_time(function() use ($matrix1) {
    ($matrix1->shape());
}, "Finding Shape Time");

//Perform subtraction first and then addition
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


measure_time(function() {

    $left = new Matrix([[2,3],[4,5]]);
    $right = new Matrix([[1,1]]);

    $left->add($right)->display();

}, "Addition Time with broadcasting");

measure_time(function() use ($matrix1, $matrix2) {
    $matrix1->sub($matrix2);
}, "Subtraction Time after Addition");

// Measure constant addition time
measure_time(function() use ($matrix1) {
    (new Matrix([[1,2,3]]))->add(0.2)->display();
}, "Addition Constant Time");

// // Measure dot product time
measure_time(function() use ($matrix1, $matrix2) {
    $matrix1->dot($matrix2);
    //var_dump($matrix1->display());
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
    $matrix1->argmax();
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
    var_dump(exp(-1));
    ((new Matrix([[-1,-2,-3,-4]]))->exp())->display();
}, "Exp Time");

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

// Measure eigen
measure_time(function() use ($matrix1,$rows,$cols) {
    $matrix1 = new Matrix([[1,2,3,4,5,6]]);
    $matrix1->reshape([6])->display();
    // $matrix1->display();
}, "reshape Time");



measure_time(function() use ($matrix1) {
    var_dump((Matrix::zeros(10,10))->shape());
}, "Zeros");

measure_time(function() use ($matrix1) {
    var_dump((Matrix::zeros_like($matrix1))->shape());
}, "Zeros like");


// Measure eigen
measure_time(function() use ($matrix1) {
    var_dump(Matrix::zeros_like($matrix1)->shape());
}, "Zeros like from passed matrix");


measure_time(function() use ($matrix1,$rows,$cols) {
    $matrix1 = new Matrix([[1,2,3,4,5,6],[7,8,9,10,11,12]]);
    $matrix1->getSlice([0,1])->display();
    // $matrix1->display();
}, "split split ==========>>");

// Measure eigen
measure_time(function() use ($matrix1) {
    ($matrix1->sum());//->display();
    // add a methode to return row col, as well as a metjhode to return the whole matrix as an array
}, "sum like");


measure_time(function() use ($matrix1) {
    ($matrix1->sqrt())->shape();
}, "sqrt like");


measure_time(function() use ($matrix1) {
    var_dump($matrix1->abs()->sum()->getData()[0][0]);
}, "get data like");

measure_time(function() use ($matrix1) {
    ($matrix1->pow(2));//->display());
}, "get data like");



measure_time(function(){

    $indices = [1, 0, 2];
    $y_pred_data = [
    [0.1, 0.9, 0.0],
    [0.2, 0.8, 0.0],
    [0.7, 0.2, 0.1]
];

    $matrix1 = new Matrix($y_pred_data);
    $result = $matrix1->getValuesFromIndices($indices);
    var_dump($result->getData()); // Output: Array ( [0] => 20 [1] => 30 [2] => nan )

}, "getValuesFromIndices");


measure_time(function() {

    // $matrix = new Matrix();
    $eye = Matrix::eye(3);

    echo "\n\n Eye Test (3X3)\n\n";

    //$eye->display();

}, "eye time");

measure_time(function() {

    $matrix = new Matrix();
    $eye = $matrix->eye(3);

    echo "\n\n Eye Test (3X3)\n\n";

    $eye->display();

}, "eye time");


measure_time(function() {

    $matrix = new Matrix();
    $matrix = $matrix->eye(3);
    echo "\n\n EYE \n\n";
    $matrix->display();
    $indices = [1, 0, 2];
    $result = $matrix->selectRowsByIndices($indices);

    echo "\n\n selectRowsByIndices \n\n";

    $result->display();

}, "eye time");




measure_time(function (){
    $matrix = new Matrix([[1, 2], [3, 4]]);
$mean_all = $matrix->mean(); // Returns 2.5
var_dump($mean_all->getData());
$mean_rows = $matrix->mean(0); // Returns [2.0, 3.0]
$mean_cols = $matrix->mean(-1); // Returns [1.5, 3.5]

echo "\n MATRIX \n";
//$matrix->display();
echo "\n MATRIX clip(0.01,0.1)\n";
$matrix->clip(0.01,0.1);///->display();

},"mean time time");


measure_time(function (){

$matrix = new Matrix([[0, -2], [3, -4]]);
$matrix1 = new Matrix([[1, -2.5], [3, -4]]);


$matrix->abs()->sub($matrix1->abs())->abs()->mul(-1)->relu(0,1)->abs();//->display();

},"sub time time");



measure_time(function() use ($matrix1) {

    $matrix1 = $matrix1->std();

    var_dump($matrix1->getData()[0][0]); // Output: Array ( [0] => 20 [1] => 30 [2] => nan )

}, "get data like");

measure_time(function (){
    $matrix = new Matrix([[1, -2], [3, 0]]);
    $where = [[1, 0], [1, 1]];
    $result = $matrix->sign();
    $result->display();
},"sign time time");

measure_time(function(){
    $matrix = new Matrix([[1, 0 ,0, 1]]);

    var_dump($matrix->sum()->getData()[0][0]);
}, "Sum along columns Time");




measure_time(function() use ($matrix1,$rows,$cols) {

    $y_true_data = [
    [1],
    [0],
    [2]
];
// $y_true->setData($y_true_data);
$y_true = new Matrix($y_true_data);

    $y_true->reshape([3]);//->display();
}, "reshape Time");



measure_time(function(){
    $matrix1 = new Matrix([[0, -2], [3, -4]]);
    $matrix2 = new Matrix([[1], [3]]);
    $matrix1->sub($matrix2);//->display();
}, "Subtraction Time");

measure_time(function(){
    $matrix1 = new Matrix([[0.5, -2], [3, -4]]);
    $matrix2 = new Matrix([[1], [3]]);
    $matrix1->relu(0.5,-100);//->display();
}, "ones_like Time");

measure_time(function(){
    $matrix1 = new Matrix([[0.5], [3]]);
    //$matrix2 = new Matrix([[1], [3]]);
    Matrix::ones_like($matrix1)->display();
}, "ones_like Time");
// measure_time(function() use ($matrix1) {

//     $matrix1->copy();//->display();
    
// }, "Copy Time");



// Test getItem
measure_time(function() use ($matrix1) {
    $originalValue = $matrix1->getItem(0, 0);
    echo "Original value at (0,0): $originalValue\n";
}, "getItem Time");

// Test setItem
measure_time(function() use ($matrix1) {
    $newValue = 999.999;
    $matrix1->setItem(0, 0, $newValue);
    echo "Set value at (0,0) to: $newValue\n";
}, "setItem Time");

// Test setItem
measure_time(function() use ($matrix1) {
    $newValue = 999.999;
    $matrix1->setItem(0, 0, $newValue);
    $modifiedValue = $matrix1->getItem(0, 0);
echo "Modified value at (0,0): $modifiedValue == $newValue\n";
}, "setItem & GET Time");


// Test setItem
measure_time(function() {
$matrix = new Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
$row = $matrix->getRow(0);  // Gets the second row (index 1)
$row->display();  // Should output: 4 5 6
}, "GET ROW Time");

measure_time(function() {
    $matrix = new Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
    $col = $matrix->getCol(1);  // Gets the second column (index 1)
    echo "Column 1:\n";
    $col->display();  // Should output: 2 5 8
}, "GET COL Time");

// Test setCol
measure_time(function() {
    $matrix = new Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
    $newCol = new Matrix([[10], [11], [12]]);
    $matrix->setCol(1, $newCol);
    echo "Matrix after setting column 1:\n";
    $matrix->display();  // Should output: 1 10 3 \n 4 11 6 \n 7 12 9
}, "SET COL Time");

// Test getRow (for comparison)
measure_time(function() {
    $matrix = new Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
    $row = $matrix->getRow(1);  // Gets the second row (index 1)
    echo "Row 1:\n";
    $row->display();  // Should output: 4 5 6
}, "GET ROW Time");

// Test setRow (for comparison)
measure_time(function() {
    $matrix = new Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
    $newRow = new Matrix([[10, 11, 12]]);
    $matrix->setRow(1, $newRow);
    echo "Matrix after setting row 1:\n";
    $matrix->display();  // Should output: 1 2 3 \n 10 11 12 \n 7 8 9
}, "SET ROW Time");



// Assuming you have a matrix with 1 row and 5 columns
$original_matrix = new Matrix([[1, 2, 3, 4, 5]]);

// Reshape it to 5 rows and 1 column
$reshaped_matrix = $original_matrix->reshape([5, 1]);

// Display the result
$reshaped_matrix->display();


// Test glorot_uniform initialization
measure_time(function() {
    $fan_in = 1000;
    $fan_out = 1000;
    
    echo "Testing glorot_uniform initialization ($fan_in, $fan_out):\n";
    
    $matrix = new Matrix();
    $glorot_matrix = Matrix::glorot_uniform($fan_in, $fan_out);
    
    // Check dimensions
    $shape = $glorot_matrix->shape();
    echo "Shape: " . implode(" x ", $shape) . "\n";
    
    // Check mean (should be close to 0)
    $mean = $glorot_matrix->mean()->getData()[0][0];
    echo "Mean: $mean\n";
    
    // Check standard deviation (should be close to sqrt(2 / (fan_in + fan_out)))
    $std = $glorot_matrix->std()->getData()[0][0];
    $expected_std = sqrt(2 / ($fan_in + $fan_out));
    echo "Standard Deviation: $std (Expected: $expected_std)\n";
    
    // Check min and max values
    // $min = $glorot_matrix->min()->getData()[0][0];
    // $max = $glorot_matrix->max()->getData()[0][0];
    // $expected_limit = sqrt(6 / ($fan_in + $fan_out));
    // echo "Min: $min, Max: $max (Expected range: -$expected_limit to $expected_limit)\n";
    
    // Optionally, you can display a small subset of the matrix
    // echo "Sample of the matrix (5x5):\n";
    // $glorot_matrix->getItem(0, 0, 5, 5)->display();
    
    return $glorot_matrix;
}, "Glorot Uniform Initialization Time");


measure_time(function() {
    $matrix = new Matrix([[1, 2], [3, 4]]);
    
    echo "Global minimum:\n";
    $matrix->min(-1, PHP_FLOAT_MAX, null)->display();

    echo "Column-wise minimum:\n";
    $matrix->min(0, PHP_FLOAT_MAX, null)->display();

    echo "Row-wise minimum:\n";
    $matrix->min(1, PHP_FLOAT_MAX, null)->display();

    echo "Minimum with initial value:\n";
    $matrix->min(-1, 2, null)->display();

    echo "Minimum with where condition:\n";
    $where = new Matrix([[1, 0], [0, 1]]);
    $matrix->min(-1, PHP_FLOAT_MAX, $where)->display();

}, "GET MIN");

measure_time(function() {
    $matrix = new Matrix([[1, 2], [3, 4]]);
    
    echo "Global maximum:\n";
    $matrix->max(-1, PHP_FLOAT_MIN, null)->display();

    echo "Column-wise maximum:\n";
    $matrix->max(0, PHP_FLOAT_MIN, null)->display();

    echo "Row-wise maximum:\n";
    $matrix->max(1, PHP_FLOAT_MIN, null)->display();

    echo "Maximum with initial value:\n";
    $matrix->max(-1, 5, null)->display();

    echo "Maximum with where condition:\n";
    $where = new Matrix([[1, 0], [0, 1]]);
    $matrix->max(-1, PHP_FLOAT_MIN, $where)->display();

}, "GET MAX");


measure_time(function() {
    echo "Testing Matrix Random Function:\n\n";

    $matrix = new Matrix();

    // Test 1: Generate a 2x3 matrix with random values between 0 and 1
    echo "1. Generate a 2x3 matrix with random values between 0 and 1:\n";
    $matrix1 = Matrix::random(2, 3);
    $matrix1->display();
    echo "\n";

    // Test 2: Generate a 3x3 matrix with random values between -1 and 1
    echo "2. Generate a 3x3 matrix with random values between -1 and 1:\n";
    $matrix2 = Matrix::random(3, 3, -1, 1);
    $matrix2->display();
    echo "\n";

    // Test 3: Generate a 1x5 matrix (vector) with random integers between 1 and 10
    echo "3. Generate a 1x5 matrix (vector) with random integers between 0 and 1 with round:\n";
    $matrix3 = Matrix::random(1, 5, 0, 1);
    $matrix3->round()->display();
    echo "\n";

    // Test 4: Generate a 4x4 matrix with random values and calculate its sum
    echo "4. Generate a 4x4 matrix with random values and calculate its sum:\n";
    $matrix4 = Matrix::random(4, 4);
    $matrix4->display();
    echo "Sum of all elements: " . $matrix4->sum(-1, 0, null)->getData()[0][0] . "\n\n";

    // Test 5: Generate a 3x3 matrix with random values and find its maximum value
    echo "5. Generate a 3x3 matrix with random values and find its maximum value:\n";
    $matrix5 = Matrix::random(3, 3);
    $matrix5->display();
    echo "Maximum value: " . $matrix5->max(-1, PHP_FLOAT_MIN, null)->getData()[0][0] . "\n\n";

    // Test 6: Generate two random matrices and multiply them
    echo "6. Generate two random matrices and multiply them:\n";
    $matrixA = Matrix::random(2, 3);
    $matrixB = Matrix::random(3, 2);
    echo "Matrix A:\n";
    $matrixA->display();
    echo "Matrix B:\n";
    $matrixB->display();
    echo "Result of A * B:\n";
    $matrixA->dot($matrixB)->display();

}, "Matrix Random Function Tests");

?>