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
        this.initialized = false;
    }

    async loadMappings() {
        if (this.initialized) return true;
        
        try {
            console.log('Starting to load mapping files...');
            const [monarchResponse, omimResponse] = await Promise.all([
                fetch('/static/data/monarch-omim.json'),
                fetch('/static/data/omim-monarch.json')
            ]);
            
            if (!monarchResponse.ok || !omimResponse.ok) {
                throw new Error('One or more responses not OK');
            }

            this.monarchMappings = await monarchResponse.json();
            this.omimMappings = await omimResponse.json();
            this.initialized = true;
            
            console.log('Mappings loaded:', {
                omimCount: Object.keys(this.omimMappings).length,
                sampleOmim: Object.keys(this.omimMappings)[0]
            });
            
            return true;
        } catch (error) {
            console.error('Error loading mapping files:', error);
            this.monarchMappings = {};
            this.omimMappings = {};
            return false;
        }
    }

    getOmimToMonarch(omimId) {
        if (!this.initialized) {
            console.error('Mapper not initialized');
            return null;
        }
        
        if (!omimId) return null;
        const cleanOmimId = omimId.replace('OMIM:', '').trim();
        
        const result = this.omimMappings[cleanOmimId];
        console.log('OMIM lookup:', {
            input: omimId,
            cleaned: cleanOmimId,
            found: !!result
        });
        
        return result || null;
    }
}

//module.exports = MonarchOmimMapper;
window.MonarchOmimMapper = MonarchOmimMapper;