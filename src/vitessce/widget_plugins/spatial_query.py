import colorsys
import json
from oxc_py import transform
from ..widget import VitesscePlugin


def _build_plugin_esm(cell_type_list):
    ct_list_js = json.dumps(cell_type_list)
    js_source = f"""
function createPlugins(utilsForPlugins) {{
    const {{
        React,
        PluginFileType,
        PluginViewType,
        PluginCoordinationType,
        PluginJointFileType,
        z,
        useCoordination,
    }} = utilsForPlugins;

    const CELL_TYPE_LIST = {ct_list_js};

    function SpatialQueryView(props) {{
        const {{ coordinationScopes }} = props;
        const [{{
            queryParams,
        }}, {{
            setQueryParams,
        }}] = useCoordination(['queryParams', 'obsType'], coordinationScopes);

        const [uuid, setUuid] = React.useState(1);
        const [queryType, setQueryType] = React.useState('grid');
        const [anchorCellType, setAnchorCellType] = React.useState(CELL_TYPE_LIST[0] ?? '');
        const [maxDist, setMaxDist] = React.useState(10);
        const [minSupport, setMinSupport] = React.useState(0.5);
        const [k, setK] = React.useState(20);
        const [nPoints, setNPoints] = React.useState(1000);

        const isAnchorMode = queryType === 'anchor-type-knn' || queryType === 'anchor-type-dist';

        const onQueryTypeChange = React.useCallback((e) => {{
            const newType = e.target.value;
            setQueryType(newType);
            // Update maxDist default when switching modes
            if (newType === 'anchor-type-knn') {{
                setMaxDist(20);
            }} else {{
                setMaxDist(10);
            }}
        }}, []);

        const radiusLabel = queryType === 'anchor-type-knn' ? 'Max. Dist.' : 'Radius';

        return (
        <div className="spatial-query">
            <p>SpatialQuery Manager</p>
            <label>
                Query type&nbsp;
                <select value={{queryType}} onChange={{onQueryTypeChange}}>
                    <option value="grid">Grid-based</option>
                    <option value="rand">Random-based</option>
                    <option value="anchor-type-knn">Anchor cell - kNN</option>
                    <option value="anchor-type-dist">Anchor cell - Radius</option>
                </select>
            </label>
            <br/>
            {{isAnchorMode && (
            <label>
                Anchor cell&nbsp;
                <select value={{anchorCellType}} onChange={{e => setAnchorCellType(e.target.value)}}>
                    {{CELL_TYPE_LIST.map(ct => (
                        <option key={{ct}} value={{ct}}>{{ct}}</option>
                    ))}}
                </select>
            </label>
            )}}
            {{isAnchorMode && <br/>}}
            <label>
                {{radiusLabel}}
                <input type="range" value={{maxDist}} onChange={{e => setMaxDist(parseFloat(e.target.value))}} min={{5}} max={{30}} step={{1}} />
                {{maxDist}}
            </label>
            <br/>
            {{queryType === 'anchor-type-knn' && (
            <label>
                k (neighbors)
                <input type="range" value={{k}} onChange={{e => setK(parseFloat(e.target.value))}} min={{5}} max={{50}} step={{1}} />
                {{k}}
                <br/>
            </label>
            )}}
            {{queryType === 'rand' && (
            <label>
                Sample points
                <input type="range" value={{nPoints}} onChange={{e => setNPoints(parseFloat(e.target.value))}} min={{100}} max={{5000}} step={{100}} />
                {{nPoints}}
                <br/>
            </label>
            )}}
            <label>
                Min. Support
                <input type="range" value={{minSupport}} onChange={{e => setMinSupport(parseFloat(e.target.value))}} min={{0}} max={{1}} step={{0.01}} />
                {{minSupport}}
            </label>
            <br/>
            <button onClick={{(e) => {{
                setQueryParams({{
                    cellTypeOfInterest: isAnchorMode ? anchorCellType : null,
                    queryType,
                    maxDist,
                    minSupport,
                    k,
                    nPoints,
                    uuid,
                }});
                setUuid(uuid+1);
            }}}}>Find patterns</button>
        </div>
        );
    }}

    const pluginCoordinationTypes = [
        new PluginCoordinationType('queryParams', null, z.object({{
            cellTypeOfInterest: z.string().nullable(),
            queryType: z.enum(['grid', 'rand', 'anchor-type-knn', 'anchor-type-dist']),
            maxDist: z.number(),
            minSupport: z.number(),
            k: z.number(),
            nPoints: z.number(),
            uuid: z.number(),
        }}).partial().nullable()),
    ];

    const pluginViewTypes = [
        new PluginViewType('spatialQuery', SpatialQueryView, ['queryParams', 'obsType']),
    ];
    return {{ pluginViewTypes, pluginCoordinationTypes }};
}}
export default {{ createPlugins }};
"""
    return transform(js_source)


class SpatialQueryPlugin(VitesscePlugin):
    """
    Spatial-Query plugin view renders controls to change parameters passed to the Spatial-Query methods.
    """
    commands = {}

    def __init__(self,
                 adata,
                 spatial_key="X_spatial",
                 label_key="cell_type",
                 feature_name="gene_name",
                 if_lognorm=True,
                 ):
        """
        Construct a new Vitessce widget.

        :param adata: AnnData.
        :type adata: anndata.AnnData
        :param str spatial_key: The key in adata.obsm that contains the (x, y) coordinates of each cell. By default, "X_spatial".
        :param str label_key: The column in adata.obs that contains the cell type labels. By default, "cell_type".
        :param str feature_name: The key in adata.var that contains the gene names. By default, "gene_name".
        :param bool if_lognorm: Whether the data in adata.X need to be log-normalized. If input is count data, set to True. By default, True.

        .. code-block:: python

            from vitessce.widget_plugins import SpatialQueryPlugin

            plugin = SpatialQueryPlugin(adata, spatial_key="X_spatial", label_key="cell_type")
            # ...
            vc.widget(plugins=[plugin], remount_on_uid_change=False)
        """
        from SpatialQuery import spatial_query
        import matplotlib.pyplot as plt  # Add as dependency / optional dependency?

        self.adata = adata
        self.spatial_key = spatial_key
        self.label_key = label_key

        self.tt = spatial_query(
            adata=adata,
            dataset='test',
            spatial_key=spatial_key,
            label_key=label_key,
            feature_name=feature_name,
            leaf_size=10,
            build_gene_index=False,
            if_lognorm=if_lognorm,
            )

        self.tab20_rgb = [[int(r * 255), int(g * 255), int(b * 255)] for (r, g, b, a) in [plt.cm.tab20(i) for i in range(20)]]

        cell_type_list = adata.obs[label_key].unique().tolist()
        self.cell_type_list = cell_type_list
        self.initial_query_params = {}

        # Build ESM with cell type list embedded directly as a JS constant
        self.plugin_esm = _build_plugin_esm(cell_type_list)

        self.additional_obs_sets = {
            "version": "0.1.3",
            "tree": [
                {
                    "name": "SpatialQuery Results",
                    "children": []
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
                "path": ["SpatialQuery Results"],
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
                    "name": f"SpatialQuery Results {sq_id}",
                    "children": [

                    ]
                }
            ]
        }

        obs_set_color = []
        n_motifs = len(fp_tree)

        for motif_i, (row_i, row) in enumerate(fp_tree.iterrows()):
            try:
                motif = row["itemsets"]
            except KeyError:
                motif = row["motifs"]
            # anchor-type queries: use neighbor_id for motif cells, grid/rand use cell_id
            if "neighbor_id" in row.index:
                cell_i = row["neighbor_id"]
            else:
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

            # Assign each motif a unique color by evenly spacing hues around the color wheel
            hue = motif_i / max(n_motifs, 1)
            r, g, b = colorsys.hls_to_rgb(hue, 0.55, 0.75)
            motif_color = [int(r * 255), int(g * 255), int(b * 255)]
            obs_set_color.append({
                "color": motif_color,
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

        max_dist = query_params.get("maxDist", 10)
        min_support = query_params.get("minSupport", 0.5)
        k = query_params.get("k", 20)
        n_points = query_params.get("nPoints", 1000)
        query_type = query_params.get("queryType", "grid")
        cell_type_of_interest = query_params.get("cellTypeOfInterest", None)

        query_uuid = query_params["uuid"]

        print(query_params)

        if query_type == "rand":
            fp_tree = self.tt.find_patterns_rand(
                max_dist=max_dist,
                n_points=int(n_points),
                min_support=min_support,
                if_display=True,
                figsize=(9, 6),
                return_cellID=True,
            )
        elif query_type == "grid":
            fp_tree, grid_pos = self.tt.find_patterns_grid(
                max_dist=max_dist,
                min_support=min_support,
                if_display=True,
                figsize=(9, 6),
                return_cellID=True,
                return_grid=True,
            )
        elif query_type == "anchor-type-knn":
            fp_tree = self.tt.motif_enrichment_knn(
                ct=cell_type_of_interest,
                k=int(k),
                max_dist=max_dist,
                min_support=min_support,
                return_cellID=True,
            )
        elif query_type == "anchor-type-dist":
            fp_tree = self.tt.motif_enrichment_dist(
                ct=cell_type_of_interest,
                max_dist=max_dist,
                min_support=min_support,
                return_cellID=True,
            )

        # TODO: implement query types that are dependent on motif selection.

        # Previous values
        additional_obs_sets = prev_config["coordinationSpace"]["additionalObsSets"]["A"]
        obs_set_color = prev_config["coordinationSpace"]["obsSetColor"]["A"]

        # Perform query
        (new_additional_obs_sets, new_obs_set_color) = self.fp_tree_to_obs_sets_tree(fp_tree, query_uuid)

        new_sq_node = new_additional_obs_sets["tree"][0]
        sq_idx = next((i for i, n in enumerate(additional_obs_sets["tree"]) if n["name"].startswith("SpatialQuery Results")), None)
        if sq_idx is not None:
            additional_obs_sets["tree"][sq_idx] = new_sq_node
        else:
            additional_obs_sets["tree"].append(new_sq_node)
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
