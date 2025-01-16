// Minimal OMIM to MONDO converter helper
document.addEventListener('DOMContentLoaded', async function() {
    console.log('Initializing OMIM converter...');
    
    // Initialize mapper and wait for mappings to load
    const mapper = new MonarchOmimMapper();
    const loaded = await mapper.loadMappings();
    
    if (!loaded) {
        console.error('Failed to load mappings');
        return;
    }
    
    // Get UI elements
    const convertButton = document.getElementById('convert-omim-button');
    const omimInput = document.getElementById('omim_id');
    const mondoInput = document.getElementById('disease_URI');
    const alertBox = document.getElementById('omim-alert');
    
    // Show conversion result message
    function showMessage(message, isSuccess = false) {
        alertBox.textContent = message;
        alertBox.className = `alert alert-${isSuccess ? 'success' : 'danger'} mt-2`;
        alertBox.style.display = 'block';
    }
    
    // Handle conversion button click
    convertButton.addEventListener('click', function() {
        console.log('Convert button clicked');
        const omimId = omimInput.value.trim();
        
        if (!omimId) {
            showMessage('Please enter an OMIM ID');
            return;
        }
        
        if (!omimId.match(/^\d+$/)) {
            showMessage('Invalid OMIM ID. Please enter numbers only');
            return;
        }
        
        console.log('Converting OMIM:', omimId);
        const result = mapper.getOmimToMonarch(omimId);
        
        if (result) {
            console.log('Conversion successful:', result);
            mondoInput.value = result.monarchId;
            showMessage(`Match for '${result.name}' found`, true);
        } else {
            console.log('No conversion found');
            mondoInput.value = '';
            showMessage('No matching MONDO ID found');
}
    });
});