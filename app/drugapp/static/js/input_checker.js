// custom command to check the formatting of 'Disease URI' and 'ML seed' inputs according to forms.py
document.addEventListener('DOMContentLoaded', function() {
    const form = document.querySelector('form');
    const diseaseURIInput = document.getElementById('disease_URI');
    const mlSeedInput = document.getElementById('ML_seed');
    const diseaseURIAlert = document.getElementById('mondo-alert');
    const mlSeedAlert = document.getElementById('ml-seed-alert');
    const submitButton = document.querySelector('input[type="submit"]');

    function validateDiseaseURI() {
        const diseaseURIPattern = /^MONDO:\d{7}$/;
        const isValid = diseaseURIPattern.test(diseaseURIInput.value);
        
        if (!isValid) {
            diseaseURIInput.setCustomValidity('Disease URI must be in format "MONDO:" followed by exactly 7 digits');
            diseaseURIInput.reportValidity();
            
            // Show alert
            diseaseURIAlert.textContent = 'Disease URI must be in format "MONDO:" followed by exactly 7 digits';
            diseaseURIAlert.className = 'alert alert-danger mt-2';
            diseaseURIAlert.style.display = 'block';
            
            return false;
        } else {
            diseaseURIInput.setCustomValidity('');
            
            // Hide alert
            diseaseURIAlert.textContent = '';
            diseaseURIAlert.className = 'alert mt-2';
            diseaseURIAlert.style.display = 'none';
            
            return true;
        }
    }

    function validateMLSeed() {
        const mlSeedValue = mlSeedInput.value.trim();
        
        if (mlSeedValue.toLowerCase() === 'random') {
            mlSeedInput.setCustomValidity('');
            
            // Hide alert
            mlSeedAlert.textContent = '';
            mlSeedAlert.className = 'alert mt-2';
            mlSeedAlert.style.display = 'none';
            
            return true;
        }
        
        // Check if it's a valid integer
        const integerPattern = /^-?\d+$/;
        const isValid = integerPattern.test(mlSeedValue);
        
        if (!isValid) {
            mlSeedInput.setCustomValidity('ML Seed must be either "random" or an integer');
            mlSeedInput.reportValidity();
            
            // Show alert
            mlSeedAlert.textContent = 'ML Seed must be either "random" or an integer';
            mlSeedAlert.className = 'alert alert-danger mt-2';
            mlSeedAlert.style.display = 'block';
            
            return false;
        } else {
            mlSeedInput.setCustomValidity('');
            
            // Hide alert
            mlSeedAlert.textContent = '';
            mlSeedAlert.className = 'alert mt-2';
            mlSeedAlert.style.display = 'none';
            
            return true;
        }
    }

    // Add real-time validation on input
    diseaseURIInput.addEventListener('input', validateDiseaseURI);
    mlSeedInput.addEventListener('input', validateMLSeed);

    // Validation on form submission
    form.addEventListener('submit', function(event) {
        const isURIValid = validateDiseaseURI();
        const isSeedValid = validateMLSeed();

        if (!isURIValid || !isSeedValid) {
            event.preventDefault(); // prevents form submission if validation fails
        }
    });
});