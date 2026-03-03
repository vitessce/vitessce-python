import base64
import colorsys
import io
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
        invokeCommand,
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
            {{queryType === 'anchor-type-knn' && (
            <label>
                k (neighbors)
                <input type="range" value={{k}} onChange={{e => setK(parseFloat(e.target.value))}} min={{5}} max={{100}} step={{1}} />
                {{k}}
                <br/>
            </label>
            )}}
            <label>
                {{radiusLabel}}
                <input type="range" value={{maxDist}} onChange={{e => setMaxDist(parseFloat(e.target.value))}} min={{5}} max={{30}} step={{1}} />
                {{maxDist}}
            </label>
            <br/>
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

    function HeatmapView(props) {{
        const {{ coordinationScopes }} = props;
        const [{{ queryParams }}] = useCoordination(['queryParams'], coordinationScopes);

        const [imgSrc, setImgSrc] = React.useState(null);
        const [loading, setLoading] = React.useState(false);

        const uuid = queryParams?.uuid;

        React.useEffect(() => {{
            if (uuid == null) return;
            setLoading(true);
            invokeCommand('get_heatmap', {{ uuid }}, []).then(([result]) => {{
                if (result?.img) {{
                    setImgSrc('data:image/png;base64,' + result.img);
                }} else {{
                    setImgSrc(null);
                }}
                setLoading(false);
            }}).catch(() => setLoading(false));
        }}, [uuid]);

        return (
        <div className="spatial-query-heatmap" style={{{{ width: '100%', height: '100%', overflow: 'auto', display: 'flex', alignItems: 'center', justifyContent: 'center' }}}}>
            {{loading && <p>Loading heatmap...</p>}}
            {{!loading && imgSrc && <img src={{imgSrc}} style={{{{ maxWidth: '100%', maxHeight: '100%' }}}} />}}
            {{!loading && !imgSrc && <p style={{{{ color: '#888' }}}}>Run a query to see the heatmap.</p>}}
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
        new PluginViewType('spatialQueryHeatmap', HeatmapView, ['queryParams']),
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
            "tree": []
        }

        self.obs_set_color = [
            {
                "color": [255, 255, 255],
                "path": ["Cell Type"],
            },
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

        self._last_fp_tree = None
        self._last_query_type = None
        self._last_cell_type_of_interest = None

        self.commands = {"get_heatmap": self._handle_get_heatmap}

    def _render_heatmap_base64(self, fp_tree, query_type, cell_type_of_interest):
        import matplotlib
        import matplotlib.pyplot as plt
        import seaborn as sns
        _prev_backend = matplotlib.get_backend()
        matplotlib.use("Agg")

        is_enrichment = query_type in ("anchor-type-knn", "anchor-type-dist")

        if is_enrichment:
            enrich = fp_tree.copy()
            enrich["frequency"] = enrich["n_center_motif"] / enrich["n_center"]
            enrich = enrich.sort_values(by="frequency", ascending=False)
            enrich["motif_group"] = [f"motif_{i+1}" for i in range(len(enrich))]
            expanded = enrich.explode("motifs")
            heatmap_data = expanded.pivot_table(
                index="motifs", columns="motif_group", values="frequency", aggfunc="first"
            )
            # Sort columns by descending frequency of each motif
            col_order = enrich.sort_values("frequency", ascending=False)["motif_group"].tolist()
            heatmap_data = heatmap_data[[c for c in col_order if c in heatmap_data.columns]]
            cbar_label = "Frequency"
            if cell_type_of_interest:
                title = f"Distribution of enriched motifs around {cell_type_of_interest}"
            else:
                title = "Distribution of enriched motifs"
            xlabel = "Motifs"
        else:
            fp = fp_tree.copy()
            fp = fp.sort_values(by="support", ascending=False).reset_index(drop=True)
            fp["motif_group"] = [f"motif_{i+1}" for i in range(len(fp))]
            expanded = fp.explode("itemsets")
            heatmap_data = expanded.pivot_table(
                index="itemsets", columns="motif_group", values="support", aggfunc="first"
            )
            # Sort columns by descending support of each motif
            col_order = fp.sort_values("support", ascending=False)["motif_group"].tolist()
            heatmap_data = heatmap_data[[c for c in col_order if c in heatmap_data.columns]]
            cbar_label = "Support"
            title = "Distribution of frequent patterns"
            xlabel = "Patterns"

        fig, ax = plt.subplots(figsize=(7, 5))
        sns.heatmap(
            heatmap_data,
            cmap="GnBu",
            linewidths=0.1,
            linecolor="lightgrey",
            annot=False,
            cbar_kws={"label": cbar_label},
            ax=ax,
        )
        ax.set_title(title, fontsize=13, pad=12)
        ax.set_ylabel("")
        ax.set_xlabel(xlabel, fontsize=11)
        ax.tick_params(axis="x", rotation=90, labelsize=10)
        ax.tick_params(axis="y", rotation=0, labelsize=10)
        fig.tight_layout()

        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        matplotlib.use(_prev_backend)
        buf.seek(0)
        return base64.b64encode(buf.read()).decode("utf-8")

    def _handle_get_heatmap(self, _message, _buffers):
        if self._last_fp_tree is None:
            return {"img": None}, []
        img_b64 = self._render_heatmap_base64(
            self._last_fp_tree,
            self._last_query_type,
            self._last_cell_type_of_interest,
        )
        return {"img": img_b64}, []

    def get_matching_cell_ids(self, cell_type, cell_i):
        cell_ids = [self.cell_i_to_cell_id[i] for i in cell_i]
        matches = []
        for cell_id in cell_ids:
            cell_ct = self.cell_id_to_cell_type[cell_id]
            if cell_ct == cell_type:
                matches.append([cell_id, None])
        return matches

    def fp_tree_to_obs_sets_tree(self, fp_tree, sq_id, anchor_ct=None):
        sq_motif_name = f"SpatialQuery Results {sq_id} — By Motif"
        sq_ct_name = f"SpatialQuery Results {sq_id} — By Cell Type"

        # Pass 1: collect per-motif data and accumulate deduped cell ids per cell type
        motif_rows = []
        ct_to_cell_ids = {}

        for _, row in fp_tree.iterrows():
            try:
                motif = row["itemsets"]
            except KeyError:
                motif = row["motifs"]
            cell_i = row["neighbor_id"] if "neighbor_id" in row.index else row["cell_id"]
            center_i = row["center_id"] if "center_id" in row.index else None
            motif_rows.append((motif, cell_i, center_i))
            for cell_type in motif:
                matching = {i for i in cell_i if self.cell_id_to_cell_type.get(self.cell_i_to_cell_id.get(i)) == cell_type}
                ct_to_cell_ids.setdefault(cell_type, set()).update(matching)

        n_motifs = len(motif_rows)
        obs_set_color = []

        # Node 1: "By Motif" — each motif is a leaf, colored by motif index
        by_motif_children = []
        for motif_i, (motif, cell_i, _center_i) in enumerate(motif_rows):
            motif_name = str(list(motif))
            hue = motif_i / max(n_motifs, 1)
            r, g, b = colorsys.hls_to_rgb(hue, 0.55, 0.75)
            motif_color = [int(r * 255), int(g * 255), int(b * 255)]
            motif_cell_types = set(motif)
            motif_cell_ids = list({
                i for i in cell_i
                if self.cell_id_to_cell_type.get(self.cell_i_to_cell_id.get(i)) in motif_cell_types
            })
            by_motif_children.append({
                "name": motif_name,
                "set": [[self.cell_i_to_cell_id[i], None] for i in motif_cell_ids if i in self.cell_i_to_cell_id]
            })
            obs_set_color.append({"color": motif_color, "path": [sq_motif_name, motif_name]})

        # Node 2: "By Cell Type" — each cell type is a leaf (union across motifs), colored by cell type
        by_ct_children = []
        for cell_type, cell_ids_set in ct_to_cell_ids.items():
            by_ct_children.append({
                "name": cell_type,
                "set": [[self.cell_i_to_cell_id[i], None] for i in sorted(cell_ids_set) if i in self.cell_i_to_cell_id]
            })
            obs_set_color.append({"color": self.ct_to_color[cell_type], "path": [sq_ct_name, cell_type]})

        # Nodes 3+: one top-level node per motif, children are cell types within that motif.
        # For anchor queries, also include an "{anchor_ct} (anchor)" child for center cells.
        per_motif_nodes = []
        for motif_i, (motif, cell_i, center_i) in enumerate(motif_rows):
            motif_name = str(list(motif))
            node_name = f"SpatialQuery Results {sq_id} \u2014 Motif {motif_i + 1}: {motif_name}"
            motif_cell_types = set(motif)
            motif_ct_children = []
            for cell_type in motif_cell_types:
                ct_ids_in_motif = list(dict.fromkeys(
                    i for i in cell_i
                    if self.cell_id_to_cell_type.get(self.cell_i_to_cell_id.get(i)) == cell_type
                ))
                motif_ct_children.append({
                    "name": cell_type,
                    "set": [[self.cell_i_to_cell_id[i], None] for i in ct_ids_in_motif if i in self.cell_i_to_cell_id]
                })
                obs_set_color.append({
                    "color": self.ct_to_color[cell_type],
                    "path": [node_name, cell_type]
                })
            # Add anchor child node if this is an anchor-mode query
            if anchor_ct is not None and center_i is not None and len(center_i) > 0:
                anchor_label = f"{anchor_ct} (anchor)"
                anchor_ids = list(dict.fromkeys(int(i) for i in center_i if i in self.cell_i_to_cell_id))
                motif_ct_children.append({
                    "name": anchor_label,
                    "set": [[self.cell_i_to_cell_id[i], None] for i in anchor_ids]
                })
                obs_set_color.append({
                    "color": self.ct_to_color.get(anchor_ct, [200, 200, 200]),
                    "path": [node_name, anchor_label]
                })
            per_motif_nodes.append({"name": node_name, "children": motif_ct_children})
            obs_set_color.append({"color": [255, 255, 255], "path": [node_name]})

        additional_obs_sets = {
            "version": "0.1.3",
            "tree": [
                {"name": sq_motif_name, "children": by_motif_children},
                {"name": sq_ct_name, "children": by_ct_children},
            ] + per_motif_nodes
        }

        obs_set_color.insert(0, {"color": [255, 255, 255], "path": [sq_motif_name]})
        obs_set_color.insert(0, {"color": [255, 255, 255], "path": [sq_ct_name]})

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
            fp_tree = fp_tree[fp_tree["if_significant"]]
        elif query_type == "anchor-type-dist":
            fp_tree = self.tt.motif_enrichment_dist(
                ct=cell_type_of_interest,
                max_dist=max_dist,
                min_support=min_support,
                return_cellID=True,
            )
            fp_tree = fp_tree[fp_tree["if_significant"]]

        # Cache for heatmap rendering
        self._last_fp_tree = fp_tree
        self._last_query_type = query_type
        self._last_cell_type_of_interest = cell_type_of_interest

        # Previous values
        additional_obs_sets = prev_config["coordinationSpace"]["additionalObsSets"]["A"]
        obs_set_color = prev_config["coordinationSpace"]["obsSetColor"]["A"]

        # Perform query
        anchor_ct = cell_type_of_interest if query_type in ("anchor-type-knn", "anchor-type-dist") else None
        (new_additional_obs_sets, new_obs_set_color) = self.fp_tree_to_obs_sets_tree(fp_tree, query_uuid, anchor_ct=anchor_ct)

        # Replace any existing SpatialQuery Results nodes (both By Motif and By Cell Type)
        existing_tree = additional_obs_sets["tree"]
        existing_tree = [n for n in existing_tree if not n["name"].startswith("SpatialQuery Results")]
        existing_tree += new_additional_obs_sets["tree"]
        additional_obs_sets["tree"] = existing_tree
        prev_config["coordinationSpace"]["additionalObsSets"]["A"] = additional_obs_sets

        obs_set_color += new_obs_set_color
        prev_config["coordinationSpace"]["obsSetColor"]["A"] = obs_set_color

        # Default selection: all motif leaf nodes under "By Motif" node
        sq_motif_node = new_additional_obs_sets["tree"][0]  # "...By Motif"
        new_obs_set_selection = [
            [sq_motif_node["name"], motif["name"]]
            for motif in sq_motif_node["children"]
        ]
        prev_config["coordinationSpace"]["obsSetSelection"]["A"] = new_obs_set_selection

        prev_config["coordinationSpace"]["obsSetExpansion"]["A"] = [[sq_motif_node["name"]]]

        return {**prev_config, "uid": f"with_query_{query_uuid}"}

    def on_config_change(self, new_config):
        query_params = new_config["coordinationSpace"]["queryParams"]["A"]
        if query_params and "uuid" in query_params:
            print(query_params)
            query_uuid = query_params.get("uuid", None)
            if new_config["uid"] != f"with_query_{query_uuid}":
                return self.run_sq(new_config)
        return None
