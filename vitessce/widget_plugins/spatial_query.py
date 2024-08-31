from esbuild_py import transform
from ..widget import VitesscePlugin


PLUGIN_ESM = transform("""
function createPlugins(utilsForPlugins) {
    const {
        React,
        PluginFileType,
        PluginViewType,
        PluginCoordinationType,
        PluginJointFileType,
        z,
        useCoordination,
    } = utilsForPlugins;
    function SpatialQueryView(props) {
        const { coordinationScopes } = props;
        const [{
            queryParams,
            obsSetSelection,
        }, {
            setQueryParams,
        }] = useCoordination(['queryParams', 'obsSetSelection', 'obsType'], coordinationScopes);

        const [uuid, setUuid] = React.useState(1);
        const [queryType, setQueryType] = React.useState('grid');
        const [maxDist, setMaxDist] = React.useState(100);
        const [minSize, setMinSize] = React.useState(4);
        const [minCount, setMinCount] = React.useState(10);
        const [minSupport, setMinSupport] = React.useState(0.5);

        const cellTypeOfInterest = obsSetSelection?.length === 1 && obsSetSelection[0][0] === "Cell Type"
            ? obsSetSelection[0][1]
            : null;

        const onQueryTypeChange = React.useCallback((e) => {
            setQueryType(e.target.value);
        }, []);

        return (
        <div className="spatial-query">
            <p>Spatial Query Manager</p>
            <label>
                Query type&nbsp;
                <select onChange={onQueryTypeChange}>
                    <option value="grid">Grid-based</option>
                    <option value="rand">Random-based</option>
                    <option value="ct-center" disabled={cellTypeOfInterest === null}>Cell type of interest</option>
                </select>
            </label>
            <br/>
            <label>
                {/* Maximum distance to consider a cell as a neighbor. */}
                Max. Dist.
                <input type="range" value={maxDist} onChange={e => setMaxDist(parseFloat(e.target.value))} min={50} max={250} step={1} />
                {maxDist}
            </label>
            <br/>
            <label>
                {/* Minimum neighborhood size for each point to consider. */}
                Min. Size
                <input type="range" value={minSize} onChange={e => setMinSize(parseFloat(e.target.value))} min={0} max={20} step={1} />
                {minSize}
            </label>
            <br/>
            <label>
                {/* Minimum number of cell type to consider. */}
                Min. Count
                <input type="range" value={minCount} onChange={e => setMinCount(parseFloat(e.target.value))} min={0} max={30} step={1} />
                {minCount}
            </label>
            <br/>
            <label>
                {/* Threshold of frequency to consider a pattern as a frequent pattern. */}
                Min. Support
                <input type="range" value={minSupport} onChange={e => setMinSupport(parseFloat(e.target.value))} min={0} max={1} step={0.01} />
                {minSupport}
            </label>
            <br/>
            {/* TODO: disDuplicates: Distinguish duplicates in patterns. */}
            <button onClick={(e) => {
                setQueryParams({
                    cellTypeOfInterest,
                    queryType,
                    maxDist,
                    minSize,
                    minCount,
                    minSupport,
                    uuid,
                });
                setUuid(uuid+1);
            }}>Find patterns</button>
        </div>
        );
    }

    const pluginCoordinationTypes = [
        new PluginCoordinationType('queryParams', null, z.object({
            cellTypeOfInterest: z.string().nullable(),
            queryType: z.enum(['grid', 'rand', 'ct-center']),
            maxDist: z.number(),
            minSize: z.number(),
            minCount: z.number(),
            minSupport: z.number(),
            disDuplicates: z.boolean(),
            uuid: z.number(),
        }).partial().nullable()),
    ];

    const pluginViewTypes = [
        new PluginViewType('spatialQuery', SpatialQueryView, ['queryParams', 'obsSetSelection', 'obsType']),
    ];
    return { pluginViewTypes, pluginCoordinationTypes };
}
export default { createPlugins };
""")


class SpatialQueryPlugin(VitesscePlugin):
    """
    Spatial-Query plugin view renders controls to change parameters passed to the Spatial-Query methods.
    """
    plugin_esm = PLUGIN_ESM
    commands = {}

    def __init__(self, adata, spatial_key="X_spatial", label_key="cell_type"):
        """
        Construct a new Vitessce widget.

        :param adata: AnnData.
        :type adata: anndata.AnnData
        :param str spatial_key: The key in adata.obsm that contains the (x, y) coordinates of each cell. By default, "X_spatial".
        :param str label_key: The column in adata.obs that contains the cell type labels. By default, "cell_type".

        .. code-block:: python

            from vitessce.widget_plugins import SpatialQueryPlugin

            plugin = SpatialQueryPlugin(adata, spatial_key="X_spatial", label_key="cell_type")
            # ...
            vc.widget(plugins=[plugin], remount_on_uid_change=False)
        """
        from SpatialQuery.spatial_query import spatial_query
        import matplotlib.pyplot as plt  # Add as dependency / optional dependency?

        self.adata = adata
        self.spatial_key = spatial_key
        self.label_key = label_key

        self.tt = spatial_query(adata=adata, dataset='test', spatial_key=spatial_key, label_key=label_key, leaf_size=10)

        self.tab20_rgb = [[int(r * 255), int(g * 255), int(b * 255)] for (r, g, b, a) in [plt.cm.tab20(i) for i in range(20)]]

        self.additional_obs_sets = {
            "version": "0.1.3",
            "tree": [
                {
                    "name": "Spatial-Query Results",
                    "children": [

                    ]
                }
            ]
        }

        self.obs_set_color = [
            {
                "color": [255, 255, 255],
                "path": ["Cell Type"],
            },
            {
                "color": [255, 255, 255],
                "path": ["Spatial-Query Results"],
            }
        ]

        self.ct_to_color = dict()

        for ct_i, cell_type in enumerate(adata.obs[label_key].unique().tolist()):
            color = self.tab20_rgb[ct_i % 20]
            self.ct_to_color[cell_type] = color
            path = ["Cell Type", cell_type]
            self.obs_set_color.append({
                "color": color,
                "path": path
            })

        self.cell_i_to_cell_id = dict(zip(range(adata.obs.shape[0]), adata.obs.index.tolist()))
        self.cell_id_to_cell_type = dict(zip(adata.obs.index.tolist(), adata.obs[label_key].tolist()))

    def get_matching_cell_ids(self, cell_type, cell_i):
        cell_ids = [self.cell_i_to_cell_id[i] for i in cell_i]
        matches = []
        for cell_id in cell_ids:
            cell_ct = self.cell_id_to_cell_type[cell_id]
            if cell_ct == cell_type:
                matches.append([cell_id, None])
        return matches

    def fp_tree_to_obs_sets_tree(self, fp_tree, sq_id):
        additional_obs_sets = {
            "version": "0.1.3",
            "tree": [
                {
                    "name": f"Spatial-Query Results {sq_id}",
                    "children": [

                    ]
                }
            ]
        }

        obs_set_color = []

        for row_i, row in fp_tree.iterrows():
            try:
                motif = row["itemsets"]
            except KeyError:
                motif = row["motifs"]
            cell_i = row["cell_id"]

            motif_name = str(list(motif))

            additional_obs_sets["tree"][0]["children"].append({
                "name": motif_name,
                "children": [
                    {
                        "name": cell_type,
                        "set": self.get_matching_cell_ids(cell_type, cell_i)
                    }
                    for cell_type in motif
                ]
            })

            obs_set_color.append({
                "color": [255, 255, 255],
                "path": [additional_obs_sets["tree"][0]["name"], motif_name]
            })

            for cell_type in motif:
                color = self.ct_to_color[cell_type]
                path = [additional_obs_sets["tree"][0]["name"], motif_name, cell_type]
                obs_set_color.append({
                    "color": color,
                    "path": path
                })
        return (additional_obs_sets, obs_set_color)

    def run_sq(self, prev_config):
        query_params = prev_config["coordinationSpace"]["queryParams"]["A"]

        max_dist = query_params.get("maxDist", 150)
        min_size = query_params.get("minSize", 4)
        # min_count = query_params.get("minCount", 10)
        min_support = query_params.get("minSupport", 0.5)
        # dis_duplicates = query_params.get("disDuplicates", False)  # if distinguish duplicates of cell types in neighborhood
        query_type = query_params.get("queryType", "grid")
        cell_type_of_interest = query_params.get("cellTypeOfInterest", None)

        query_uuid = query_params["uuid"]

        params_dict = dict(
            max_dist=max_dist,
            min_size=min_size,
            # min_count=min_count,
            min_support=min_support,
            # dis_duplicates=dis_duplicates,
            if_display=True,
            fig_size=(9, 6),
            return_cellID=True,
        )
        print(params_dict)

        # TODO: add unit tests for this functionality

        if query_type == "rand":
            # TODO: implement param similar to return_grid for find_patterns_rand (to return the random points used)
            fp_tree = self.tt.find_patterns_rand(**params_dict)
        elif query_type == "grid":
            params_dict["return_grid"] = True
            fp_tree, grid_pos = self.tt.find_patterns_grid(**params_dict)
        elif query_type == "ct-center":
            fp_tree = self.tt.motif_enrichment_knn(
                ct=cell_type_of_interest,
                k=20,  # TODO: make this a parameter in the UI.
                min_support=min_support,
                # dis_duplicates=dis_duplicates,
                return_cellID=True,
            )
            print(fp_tree)

        # TODO: implement query types that are dependent on motif selection.

        # Previous values
        additional_obs_sets = prev_config["coordinationSpace"]["additionalObsSets"]["A"]
        obs_set_color = prev_config["coordinationSpace"]["obsSetColor"]["A"]

        # Perform query
        (new_additional_obs_sets, new_obs_set_color) = self.fp_tree_to_obs_sets_tree(fp_tree, query_uuid)

        additional_obs_sets["tree"][0] = new_additional_obs_sets["tree"][0]
        prev_config["coordinationSpace"]["additionalObsSets"]["A"] = additional_obs_sets

        obs_set_color += new_obs_set_color
        prev_config["coordinationSpace"]["obsSetColor"]["A"] = obs_set_color

        motif_to_select = new_additional_obs_sets["tree"][0]["children"][0]["name"]
        new_obs_set_selection = [[new_additional_obs_sets["tree"][0]["name"], motif_to_select, node["name"]] for node in new_additional_obs_sets["tree"][0]["children"][0]["children"]]
        prev_config["coordinationSpace"]["obsSetSelection"]["A"] = new_obs_set_selection

        # TODO: need to fix bug that prevents this from working
        # Reference: https://github.com/vitessce/vitessce/blob/774328ab5c4436576dd2e8e4fff0758d6c6cce89/packages/view-types/obs-sets-manager/src/ObsSetsManagerSubscriber.js#L104
        prev_config["coordinationSpace"]["obsSetExpansion"]["A"] = [path[:-1] for path in new_obs_set_selection]

        return {**prev_config, "uid": f"with_query_{query_uuid}"}

    def on_config_change(self, new_config):
        query_params = new_config["coordinationSpace"]["queryParams"]["A"]
        if query_params and "uuid" in query_params:
            print(query_params)
            query_uuid = query_params.get("uuid", None)
            if new_config["uid"] != f"with_query_{query_uuid}":
                return self.run_sq(new_config)
        return None
