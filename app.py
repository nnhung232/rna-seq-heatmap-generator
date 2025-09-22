import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import streamlit as st
import io
from matplotlib.colors import LinearSegmentedColormap

st.set_page_config(page_title="Heatmap Generator", layout="wide")

# Create two columns for layout
col1, col2 = st.columns([1, 3])  # Adjust ratio as needed


with col1:
    st.header("Heatmap Generator")
    # Upload Excel File
    uploaded_first_file = st.file_uploader("Upload your first Excel file", type=["xlsx"])
    uploaded_second_file = st.file_uploader("Upload your second Excel file", type=["xlsx"])

    # Sidebar controls
    st.sidebar.header("Controls")
    sheet_first = "sigDEG_FC1"
    sheet_second = "sigDEG_FC1"
    if uploaded_first_file is not None:
        xl_first = pd.ExcelFile(uploaded_first_file)
        sheet_first = st.sidebar.selectbox("Select First Sheet", xl_first.sheet_names, index=xl_first.sheet_names.index("sigDEG_FC1") if "sigDEG_FC1" in xl_first.sheet_names else 0)
    if uploaded_second_file is not None:
        xl_second = pd.ExcelFile(uploaded_second_file)
        sheet_second = st.sidebar.selectbox("Select Second Sheet", xl_second.sheet_names, index=xl_second.sheet_names.index("sigDEG_FC1") if "sigDEG_FC1" in xl_second.sheet_names else 0)

    # Colormap selection
    colormap_options = {
        "Blue-Black-Yellow": ["blue", "black", "yellow"],
        "Red-White-Blue": ["red", "white", "blue"],
        "Viridis": "viridis",
        "Plasma": "plasma",
        "Custom": st.sidebar.text_input("Custom Colormap (comma separated)", value="blue,black,yellow")
    }
    cmap_choice = st.sidebar.selectbox("Colormap", list(colormap_options.keys()), index=0)
    cmap_value = colormap_options[cmap_choice]
    if isinstance(cmap_value, str) and cmap_choice == "Custom":
        cmap_value = [c.strip() for c in cmap_value.split(",")]

    # vmin/vmax selection
    vmin = st.sidebar.number_input("vmin", value=None, format="%g")
    vmax = st.sidebar.number_input("vmax", value=None, format="%g")

    # Map order selection
    # map_order = st.sidebar.radio("Map Order", ["Leaf-Root", "Root-Leaf"], index=0)

    # Title, labels, font sizes
    plot_title = st.sidebar.text_input("Plot Title", value="Heatmap (All Data)")
    xlabel = st.sidebar.text_input("X Label", value="Tissues / Conditions")
    ylabel = st.sidebar.text_input("Y Label", value="")
    fontsize_title = st.sidebar.number_input("Title Fontsize", value=14, min_value=8, max_value=32)
    fontsize_xtick = st.sidebar.number_input("Xtick Fontsize", value=10, min_value=3, max_value=24)
    fontsize_ytick = st.sidebar.number_input("Ytick Fontsize", value=6, min_value=3, max_value=24)

    # Custom x tick labels
    xtick_label_1 = st.sidebar.text_input("Xtick Label 1", value="logFC_first")
    xtick_label_2 = st.sidebar.text_input("Xtick Label 2", value="logFC_second")

    if uploaded_first_file is not None and uploaded_second_file is not None:
        st.success("Files uploaded successfully!")
        if st.button("Generate Heatmap"):
            try:
                # === Step 1. Read Excel File ===
                df_first = pd.read_excel(uploaded_first_file, sheet_name=sheet_first, header=0)
                df_second = pd.read_excel(uploaded_second_file, sheet_name=sheet_second, header=0)
                # Check required columns
                required_cols = {"gene_id", "gene_symbol", "logFC"}
                if not required_cols.issubset(df_first.columns):
                    st.error(f"Missing required columns in First file. Please ensure your file has: {required_cols}")
                elif not required_cols.issubset(df_second.columns):
                    st.error(f"Missing required columns in Second file. Please ensure your file has: {required_cols}")
                else:
                    columns_to_keep_first = ['gene_id', 'gene_symbol'] + [col for col in df_first.columns if str(col).startswith('logFC')]
                    df_first = df_first[columns_to_keep_first]
                    columns_to_keep_second = ['gene_id', 'gene_symbol'] + [col for col in df_second.columns if str(col).startswith('logFC')]
                    df_second = df_second[columns_to_keep_second]
                    st.write("Data Preview:")
                    st.dataframe(df_first)
                    st.dataframe(df_second)
                    # Left join first and second dataframes on gene_id
                    df = pd.merge(df_first, df_second, on=['gene_id', 'gene_symbol'], how='left', suffixes=('_first', '_second'))
                    st.dataframe(df)
                    df["logFC_first"] = df["logFC_first"].fillna(0)
                    df["logFC_second"] = df["logFC_second"].fillna(0)
                    st.dataframe(df)
                    default_max_abs = max(abs(df["logFC_first"].min()), abs(df["logFC_first"].max()))
                    auto_vmin = -default_max_abs
                    auto_vmax = default_max_abs
                    # Use sidebar vmin/vmax if set, else auto
                    vmin_final = vmin if vmin is not None and vmin != 0 else auto_vmin
                    vmax_final = vmax if vmax is not None and vmax != 0 else auto_vmax
                    def combine_gene_info(row):
                        if pd.notna(row['gene_symbol']) and str(row['gene_symbol']).strip() != "":
                            return f"{row['gene_id']} - {row['gene_symbol']}"
                        return str(row['gene_id'])
                    df['gene_name'] = df.apply(combine_gene_info, axis=1)
                    df_all = df.copy()
                    df_all.set_index("gene_name", inplace=True)
                    # Map order
                    # if map_order == "Leaf-Root":
                    #     df_all = df_all[["logFC_leaf", "logFC_root"]]
                    # else:
                    #     df_all = df_all[["logFC_root", "logFC_leaf"]]
                    df_all = df_all[["logFC_first", "logFC_second"]]
                    df_all.columns = [xtick_label_1, xtick_label_2]
                    df_all_sorted = df_all.sort_values(xtick_label_1, ascending=False)
                    # === Step 4. Generate Heatmaps ===
                    if isinstance(cmap_value, list):
                        custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", cmap_value)
                    else:
                        custom_cmap = cmap_value
                    def save_heatmap(data, title, vmin, vmax):
                        num_genes = len(data)
                        base_height = min(10, max(6, num_genes * 0.15))
                        plt.figure(figsize=(6, base_height))
                        sns.heatmap(
                            data,
                            cmap=custom_cmap,
                            center=0,
                            vmin=vmin,
                            vmax=vmax,
                            linewidths=0,
                            cbar_kws={'label': 'Log₂ Fold Change'},
                            yticklabels=data.index
                        )
                        plt.title(title, fontsize=fontsize_title)
                        plt.xlabel(xlabel, fontsize=fontsize_xtick)
                        plt.ylabel(ylabel, fontsize=fontsize_ytick)
                        plt.yticks(rotation=0, fontsize=fontsize_ytick)
                        plt.xticks(fontsize=fontsize_xtick)
                        plt.tight_layout()
                        buf = io.BytesIO()
                        plt.savefig(buf, format="png", dpi=300, bbox_inches="tight")
                        buf.seek(0)
                        plt.close()
                        return buf
                    st.session_state["heatmap_buffer_all"] = save_heatmap(df_all_sorted, plot_title, vmin_final, vmax_final)
                    st.success("All heatmaps generated successfully!")
                    st.download_button(
                        label="Download All Data Heatmap",
                        data=st.session_state["heatmap_buffer_all"],
                        file_name="heatmap_all.png",
                        mime="image/png"
                    )

            except Exception as e:
                st.error(f"An error occurred: {e}")
                #please print stack trace
                import traceback
                st.text(traceback.format_exc())

with col2:
    # Display all five heatmaps in a grid
    heatmaps = [
        ("Heatmap (All Data)", "heatmap_buffer_all"),
    ]
    cols = st.columns(2)
    for i, (title, key) in enumerate(heatmaps):
        with cols[i % 2]:
            if key in st.session_state:
                st.header(title)
                st.image(st.session_state[key], use_container_width=True)




