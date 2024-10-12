import numpy as np

class Loss:
    def __init__(self):
        self.trainable_layers = []
        self.accumulated_sum = 0
        self.accumulated_count = 0

    def remember_trainable_layers(self, trainable_layers):
        self.trainable_layers = trainable_layers

    def regularization_loss(self):
        regularization_loss = 0
        for layer in self.trainable_layers:
            if layer.weight_regularizer_l1 > 0:
                regularization_loss += layer.weight_regularizer_l1 * np.sum(np.abs(layer.weights))
            if layer.weight_regularizer_l2 > 0:
                regularization_loss += layer.weight_regularizer_l2 * np.sum(layer.weights ** 2)
            if layer.bias_regularizer_l1 > 0:
                regularization_loss += layer.bias_regularizer_l1 * np.sum(np.abs(layer.biases))
            if layer.bias_regularizer_l2 > 0:
                regularization_loss += layer.bias_regularizer_l2 * np.sum(layer.biases ** 2)
        return regularization_loss

    def calculate(self, output, y, include_regularization=False):
        sample_losses = self.forward(output, y)
        print("Sample losses:", sample_losses)
        data_loss = np.mean(sample_losses)

        self.accumulated_sum += np.sum(sample_losses)
        self.accumulated_count += sample_losses.shape[0]

        if not include_regularization:
            return data_loss

        regularization_loss = self.regularization_loss()
        return data_loss, regularization_loss

    def new_pass(self):
        self.accumulated_sum = 0
        self.accumulated_count = 0

    def calculate_accumulated(self, include_regularization=False):
        data_loss = self.accumulated_sum / self.accumulated_count

        if not include_regularization:
            return data_loss

        regularization_loss = self.regularization_loss()
        return data_loss, regularization_loss

class Loss_CategoricalCrossentropy(Loss):
    def forward(self, y_pred, y_true):
        samples = y_pred.shape[0]
        y_pred_clipped = np.clip(y_pred, 1e-7, 1 - 1e-7)

        if len(y_true.shape) == 1:
            correct_confidences = y_pred_clipped[range(samples), y_true]
        elif len(y_true.shape) == 2:
            correct_confidences = np.sum(y_pred_clipped * y_true, axis=1)

        negative_log_likelihoods = -np.log(correct_confidences)
        return negative_log_likelihoods

    def backward(self, dvalues, y_true):
        samples = dvalues.shape[0]
        labels = dvalues.shape[1]

        print(samples)
        
        if len(y_true.shape) == 1:
            y_true = np.eye(labels)[y_true]


        self.dinputs = -y_true / dvalues   
        self.dinputs = self.dinputs / samples
       

class Loss_BinaryCrossentropy(Loss):
    def forward(self, y_pred, y_true):
        y_pred_clipped = np.clip(y_pred, 1e-7, 1 - 1e-7)

        positive_loss = y_true * np.log(y_pred_clipped)
        negative_loss = (1 - y_true) * np.log(1 - y_pred_clipped)

        sample_losses = -(positive_loss + negative_loss)
        sample_losses = np.mean(sample_losses, axis=-1)
        return sample_losses

# Backward pass
def backward ( self , dvalues , y_true ):
# Number of samples
samples = len (dvalues)
# Number of outputs in every sample
# We'll use the first sample to count them
outputs = len (dvalues[ 0 ])
# Clip data to prevent division by 0
# Clip both sides to not drag mean towards any value
clipped_dvalues = np.clip(dvalues, 1e-7 , 1 - 1e-7 )
# Calculate gradient
self.dinputs = - (y_true / clipped_dvalues -
( 1 - y_true) / ( 1 - clipped_dvalues)) / outputs
# Normalize gradient
self.dinputs = self.dinputs / samples

class Loss_MeanSquaredError(Loss):
    def forward(self, y_pred, y_true):
        sample_losses = np.mean((y_true - y_pred) ** 2, axis=-1)
        return sample_losses

    def backward(self, dvalues, y_true):
        samples = dvalues.shape[0]
        outputs = dvalues.shape[1]

        self.dinputs = -2 * (y_true - dvalues) / outputs
        self.dinputs = self.dinputs / samples

class Loss_MeanAbsoluteError(Loss):
    def forward(self, y_pred, y_true):
        sample_losses = np.mean(np.abs(y_true - y_pred), axis=-1)
        return sample_losses

    def backward(self, dvalues, y_true):
        samples = dvalues.shape[0]
        outputs = dvalues.shape[1]

        self.dinputs = np.sign(y_true - dvalues) / outputs
        self.dinputs = self.dinputs / samples

# Test example for each loss function
def test_loss_functions():
    y_pred_data = np.array([
        [0.1, 0.9, 0.0],
        [0.2, 0.8, 0.0],
        [0.7, 0.2, 0.1]
    ])
    
    y_true_data = np.array([
        [1],
        [0],
        [2]
    ])

    # Reshape y_true to be broadcastable for MAE and MSE calculation
    y_true_broadcastable = y_true_data.reshape(-1, 1)

    # Categorical Crossentropy
    loss_cce = Loss_CategoricalCrossentropy()
    cce_loss = loss_cce.calculate(y_pred_data, y_true_data.flatten())
    print("Categorical Crossentropy Loss:", cce_loss)
    loss_cce.backward(y_pred_data, y_true_data.flatten())
    print("Backward pass (CCE):", loss_cce.dinputs)

    # Binary Crossentropy
    y_pred_binary_data = np.array([
        [0.9],
        [0.2],
        [0.7]
    ])
    
    y_true_binary_data = np.array([
        [1],
        [0],
        [1]
    ])
    loss_bce = Loss_BinaryCrossentropy()
    bce_loss = loss_bce.calculate(y_pred_binary_data, y_true_binary_data)
    print("Binary Crossentropy Loss:", bce_loss)
    loss_bce.backward(y_pred_binary_data, y_true_binary_data)
    print("Backward pass (BCE):", loss_bce.dinputs)

    # Mean Squared Error
    loss_mse = Loss_MeanSquaredError()
    mse_loss = loss_mse.calculate(y_pred_data, y_true_broadcastable)
    print("Mean Squared Error Loss:", mse_loss)
    loss_mse.backward(y_pred_data, y_true_broadcastable)
    print("Backward pass (MSE):", loss_mse.dinputs)

    # Mean Absolute Error
    loss_mae = Loss_MeanAbsoluteError()
    mae_loss = loss_mae.calculate(y_pred_data, y_true_broadcastable)
    print("Mean Absolute Error Loss:", mae_loss)
    loss_mae.backward(y_pred_data, y_true_broadcastable)
    print("Backward pass (MAE):", loss_mae.dinputs)

test_loss_functions()
