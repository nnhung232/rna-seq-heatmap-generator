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
    uploaded_guided_file = st.file_uploader("Upload your Guided Excel file ", type=["xlsx"])
    uploaded_followed_file = st.file_uploader("Upload your Followed Excel file", type=["xlsx"])

    # Sidebar controls
    st.sidebar.header("Controls")
    sheet_guided = "sigDEG_FC1"
    sheet_followed = "sigDEG_FC1"
    if uploaded_guided_file is not None:
        xl_guided = pd.ExcelFile(uploaded_guided_file)
        sheet_guided = st.sidebar.selectbox("Select Guided Sheet", xl_guided.sheet_names, index=xl_guided.sheet_names.index("sigDEG_FC1") if "sigDEG_FC1" in xl_guided.sheet_names else 0)
    if uploaded_followed_file is not None:
        xl_followed = pd.ExcelFile(uploaded_followed_file)
        sheet_followed = st.sidebar.selectbox("Select Followed Sheet", xl_followed.sheet_names, index=xl_followed.sheet_names.index("sigDEG_FC1") if "sigDEG_FC1" in xl_followed.sheet_names else 0)

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
    map_order = st.sidebar.radio("Map Order", ["Guided-Followed", "Followed-Guided"], index=0)

    # Title, labels, font sizes
    plot_title = st.sidebar.text_input("Plot Title", value="Heatmap (All Data)")
    xlabel = st.sidebar.text_input("X Label", value="Tissues / Conditions")
    ylabel = st.sidebar.text_input("Y Label", value="")
    fontsize_title = st.sidebar.number_input("Title Fontsize", value=14, min_value=8, max_value=32)
    fontsize_xtick = st.sidebar.number_input("Xtick Fontsize", value=10, min_value=3, max_value=24)
    fontsize_ytick = st.sidebar.number_input("Ytick Fontsize", value=6, min_value=3, max_value=24)

    # Custom x tick labels
    xtick_label_left = st.sidebar.text_input("Xtick Left")
    xtick_label_right = st.sidebar.text_input("Xtick Right")

    # Check if at least one file is uploaded
    if uploaded_guided_file is not None:
        if uploaded_followed_file is not None:
            st.success("Both files uploaded successfully!")
        else:
            st.success("Guided file uploaded successfully! You can generate a single-column heatmap.")
        if st.button("Generate Heatmap"):
            try:
                # === Step 1. Read Excel File ===
                df_guided = pd.read_excel(uploaded_guided_file, sheet_name=sheet_guided, header=0)
                
                # Check required columns for guided file
                required_cols = {"gene_id", "gene_symbol", "logFC"}
                if not required_cols.issubset(df_guided.columns):
                    st.error(f"Missing required columns in Guided file. Please ensure your file has: {required_cols}")
                else:
                    # Process guided file
                    columns_to_keep_guided = ['gene_id', 'gene_symbol'] + [col for col in df_guided.columns if str(col).startswith('logFC')]
                    df_guided = df_guided[columns_to_keep_guided]
                    
                    # Process followed file if uploaded
                    if uploaded_followed_file is not None:
                        df_followed = pd.read_excel(uploaded_followed_file, sheet_name=sheet_followed, header=0)
                        if not required_cols.issubset(df_followed.columns):
                            st.error(f"Missing required columns in Followed file. Please ensure your file has: {required_cols}")
                        else:
                            columns_to_keep_followed = ['gene_id', 'gene_symbol'] + [col for col in df_followed.columns if str(col).startswith('logFC')]
                            df_followed = df_followed[columns_to_keep_followed]
                            
                            st.write("Data Preview:")
                            st.dataframe(df_guided)
                            st.dataframe(df_followed)
                            
                            # Left join guided and followed dataframes on gene_id
                            df = pd.merge(df_guided, df_followed, on=['gene_id', 'gene_symbol'], how='left', suffixes=('_guided', '_followed'))
                            st.dataframe(df)
                            df["logFC_guided"] = df["logFC_guided"].fillna(0)
                            df["logFC_followed"] = df["logFC_followed"].fillna(0)
                            st.dataframe(df)
                            
                            default_max_abs = max(abs(df["logFC_guided"].min()), abs(df["logFC_guided"].max()))
                            auto_vmin = -default_max_abs
                            auto_vmax = default_max_abs
                            
                            # Create gene_name column
                            def combine_gene_info(row):
                                if pd.notna(row['gene_symbol']) and str(row['gene_symbol']).strip() != "":
                                    return f"{row['gene_id']} - {row['gene_symbol']}"
                                return str(row['gene_id'])
                            df['gene_name'] = df.apply(combine_gene_info, axis=1)
                            
                            df_all = df.copy()
                            df_all.set_index("gene_name", inplace=True)
                            
                            # Map order for two columns
                            if map_order == "Guided-Followed":
                                df_all = df_all[["logFC_guided", "logFC_followed"]]
                            else:
                                df_all = df_all[["logFC_followed", "logFC_guided"]]

                            if xtick_label_left != "" and xtick_label_right != "":
                                df_all.columns = [xtick_label_left, xtick_label_right]
                    else:
                        # Single file mode - only guided file
                        st.write("Data Preview:")
                        st.dataframe(df_guided)
                        
                        df = df_guided.copy()
                        
                        default_max_abs = max(abs(df["logFC"].min()), abs(df["logFC"].max()))
                        auto_vmin = -default_max_abs
                        auto_vmax = default_max_abs
                        
                        # Create gene_name column
                        def combine_gene_info(row):
                            if pd.notna(row['gene_symbol']) and str(row['gene_symbol']).strip() != "":
                                return f"{row['gene_id']} - {row['gene_symbol']}"
                            return str(row['gene_id'])
                        df['gene_name'] = df.apply(combine_gene_info, axis=1)
                        
                        df_all = df.copy()
                        df_all.set_index("gene_name", inplace=True)
                        df_all = df_all[["logFC"]]  # Single column
                        
                        # Use custom label if provided, otherwise use default
                        if xtick_label_left != "":
                            df_all.columns = [xtick_label_left]
                    
                    # Use sidebar vmin/vmax if set, else auto
                    vmin_final = vmin if vmin is not None and vmin != 0 else auto_vmin
                    vmax_final = vmax if vmax is not None and vmax != 0 else auto_vmax
                    
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
                    st.session_state["heatmap_buffer_all"] = save_heatmap(df_all, plot_title, vmin_final, vmax_final)
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




