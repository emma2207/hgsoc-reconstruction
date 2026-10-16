def buildModalityOutputContext(dataset, isPseudobulk, ncells, modalities) {
    def modalityGroup = isPseudobulk ? 'pseudobulk' : modalities.join('_vs_')
    def experimentPath = isPseudobulk
        ? "${dataset}/pseudobulk/ncells_${ncells}"
        : "${dataset}/${modalityGroup}/ncells_null"

    [modalityGroup: modalityGroup, experimentPath: experimentPath]
}

def buildReadDepthContext(readDepthParam, dataset, isPseudobulk, modalities, ncells = null) {
    def resolvedReadDepth =
        (readDepthParam instanceof Map)
            ? (readDepthParam[dataset] ?: readDepthParam.default ?: ["default": 0])
            : ["default": readDepthParam]

    def resolvedUpper =
        (resolvedReadDepth instanceof Map)
            ? (resolvedReadDepth.upper ?: 1000)
            : 1000

    def resolvedDefault =
        (resolvedReadDepth instanceof Map)
            ? (resolvedReadDepth.default ?: 0)
            : resolvedReadDepth

    def modalityReadDepths =
        (resolvedReadDepth instanceof Map)
            ? resolvedReadDepth.findAll { k, v -> k != 'default' && k != 'upper' }
            : resolvedReadDepth

    def datasetReadDepth =
        isPseudobulk
            ? ["default": resolvedDefault]
            : modalityReadDepths

    def modalityReadDepthTag = modalities.collect { mod ->
        "${mod}_${datasetReadDepth[mod] ?: resolvedDefault}"
    }.join('_')

    def readDepthPairs = datasetReadDepth.collect { k, v -> "${k}:${v}" }.join(',')

    def pseudobulkReadDepth =
        (datasetReadDepth instanceof Map)
            ? (datasetReadDepth.default ?: 0)
            : datasetReadDepth

    def realDataNcellsPath =
        (!isPseudobulk && dataset in ['hgsoc', 'hgsoc-new'] && modalities?.size() >= 2)
            ? "real_data/${modalities[0]}_vs_${modalities[1]}/ncells_null"
            : "real_data/ncells_null"

    def modalityOutputContext = buildModalityOutputContext(dataset, isPseudobulk, ncells, modalities)

    [
        resolvedReadDepth: resolvedReadDepth,
        resolvedDefault: resolvedDefault,
        resolvedUpper: resolvedUpper,
        datasetReadDepth: datasetReadDepth,
        modalityReadDepthTag: modalityReadDepthTag,
        readDepthPairs: readDepthPairs,
        pseudobulkReadDepth: pseudobulkReadDepth,
        realDataNcellsPath: realDataNcellsPath,
        modalityGroup: modalityOutputContext.modalityGroup,
        experimentPath: modalityOutputContext.experimentPath,
    ]
}