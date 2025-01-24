/***
OMIM CONVERTER: a JavaScript conversion tool for Open Data Science URIs
Created on January 13th 2025
@author: Niccolò Bianchi [https://github.com/NCMBianchi]

'mapper.js' is a streamline ID conversion tool to convert Monarch Initiative's
custom IDs to OMIM IDs, fetching corresponding IDs and names from mapping file
'monarch-omim.json' that can be generated/updated with 'updateMapping.py'.
***/

class MonarchOmimMapper {
    constructor() {
        this.monarchMappings = {};
        this.omimMappings = {};
        this.loadMappings();
    }

    async loadMappings() {
        try {
            // Load both mapping files
            const [monarchResponse, omimResponse] = await Promise.all([
                // modified according to this specific application
                fetch('/static/data/monarch-omim.json'),
                fetch('/static/data/omim-monarch.json')
            ]);
            
            this.monarchMappings = await monarchResponse.json();
            this.omimMappings = await omimResponse.json();
        } catch (error) {
            console.error('Error loading mapping files:', error);
            this.monarchMappings = {};
            this.omimMappings = {};
        }
    }

    getMonarchToOmim(monarchId) {
        if (!monarchId) return null;
        const cleanMonarchId = monarchId.trim().toUpperCase();
        return this.monarchMappings[cleanMonarchId] || null;
    }

    getOmimToMonarch(omimId) {
        if (!omimId) return null;
        const cleanOmimId = omimId.replace('OMIM:', '').trim();
        return this.omimMappings[cleanOmimId] || null;
    }

    // Get just the name for a Monarch Initiative ID
    getDiseaseName(monarchId) {
        const data = this.getMonarchToOmim(monarchId);
        return data ? data.name : null;
    }

    // Get just the name for an OMIM ID
    getOmimName(omimId) {
        const data = this.getOmimToMonarch(omimId);
        return data ? data.name : null;
    }
}

//module.exports = MonarchOmimMapper;
window.MonarchOmimMapper = MonarchOmimMapper;