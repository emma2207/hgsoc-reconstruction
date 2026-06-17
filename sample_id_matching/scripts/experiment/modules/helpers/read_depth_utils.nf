def buildReadDepthContext(readDepthParam, dataset, isPseudobulk, modalities) {
    def resolvedReadDepth =
        (readDepthParam instanceof Map)
            ? (readDepthParam[dataset] ?: readDepthParam.default ?: ["default": 0])
            : ["default": readDepthParam]

    def resolvedDefault =
        (resolvedReadDepth instanceof Map)
            ? (resolvedReadDepth.default ?: 0)
            : resolvedReadDepth

    def datasetReadDepth =
        isPseudobulk
            ? ["default": resolvedDefault]
            : resolvedReadDepth

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

    [
        resolvedReadDepth: resolvedReadDepth,
        resolvedDefault: resolvedDefault,
        datasetReadDepth: datasetReadDepth,
        modalityReadDepthTag: modalityReadDepthTag,
        readDepthPairs: readDepthPairs,
        pseudobulkReadDepth: pseudobulkReadDepth,
        realDataNcellsPath: realDataNcellsPath,
    ]
}