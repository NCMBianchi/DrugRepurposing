document.addEventListener('DOMContentLoaded', function() {
    var minSimField = document.getElementById('inp_minimum_sim');
    var simThresholdField = document.getElementById('sim_t');
    
    minSimField.addEventListener('input', function() {
        var minSimValue = parseFloat(this.value);
        simThresholdField.min = minSimValue;
    });
});
