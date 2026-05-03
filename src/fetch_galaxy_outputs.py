#!/usr/bin/env -S pixi exec --spec bioblend -- python

import argparse
import os
import pprint
import sys

from pathlib import Path

from bioblend.galaxy import GalaxyInstance


def get_datasets_from_history(
    galaxy: GalaxyInstance, history_name: str, dataset_name_prefix: str
):
    """Return the datasets in the history 'history_name'

    galaxy: the GalaxyInstance to use
    history_name: the name of the history to look in
    """
    histories = galaxy.histories.get_histories(name=history_name)
    if len(histories) != 1:
        raise ValueError(
            f"Expected exactly one history named '{history_name}', found {len(histories)}"
        )

    history_id = histories[0]["id"]
    datasets = galaxy.histories.show_history(
        history_id,
        contents=True,
        types=["dataset"],
    )
    matched_datasets = [
        d for d in datasets if d.get("name").startswith(dataset_name_prefix)
    ]
    if len(matched_datasets) != 1:
        raise ValueError(
            f"Expected exactly one dataset with name starting with '{dataset_name_prefix}' in history "
            f"'{history_name}', found none"
        )

    return matched_datasets[0]


def get_datasets_from_collection(
    galaxy: GalaxyInstance, history_name: str, collection_name_prefix: str
):
    """Return the datasets in the list type dataset collection whose name starts with
    collection_name_prefix which is in the history 'history_name'

    galaxy: the GalaxyInstance to use
    history_name: the name of the history to look in
    collection_name_prefix: the prefix of name of the collection to look for
    """
    histories = galaxy.histories.get_histories(name=history_name)
    if len(histories) != 1:
        raise ValueError(
            f"Expected exactly one history named '{history_name}', found {len(histories)}"
        )

    history_id = histories[0]["id"]
    collections = galaxy.histories.show_history(
        history_id,
        contents=True,
        types=["dataset_collection"],
    )
    matched_collections = [
        c for c in collections if c.get("name").startswith(collection_name_prefix)
    ]
    if len(matched_collections) != 1:
        raise ValueError(
            f"Expected exactly one collection with name starting with '{collection_name_prefix}' in history "
            f"'{history_name}', found {len(matched_collections)}"
        )

    collection_id = matched_collections[0]["id"]
    collection_name = matched_collections[0]["name"]
    collection = galaxy.dataset_collections.show_dataset_collection(collection_id)
    if collection.get("collection_type") != "list":
        raise ValueError(
            f"Collection '{collection_name}' is type "
            f"'{collection.get('collection_type')}', expected 'list'"
        )

    elements = collection.get("elements", [])
    datasets = []
    for element in elements:
        # For list collections, element objects should be datasets (HDAs).
        dataset = element.get("object")
        if dataset is None:
            continue
        datasets.append(dataset)

    return datasets


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fetch enterovirus-related outputs from a Galaxy history"
    )
    parser.add_argument(
        "--galaxy_url",
        default="https://usegalaxy.eu",
        help="URL for the Galaxy server to work on",
    )
    parser.add_argument("output_dir", help="Directory to save the output datasets to")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    if not output_dir.exists():
        print(f"Output directory '{output_dir}' does not exist, creating it")
        output_dir.mkdir(parents=True, exist_ok=True)

    api_key = os.environ.get("API_KEY")
    if api_key is None:
        exit("API_KEY needed, please set this environment variable")

    gi = GalaxyInstance(url=args.galaxy_url, key=api_key)

    datasets = get_datasets_from_collection(
        galaxy=gi,
        history_name="TZ Red Eyes 4",
        collection_name_prefix="bcftools consensus on collection",
    )

    for dataset in datasets:
        dataset_id = dataset.get("id")
        name = dataset.get("name")
        for tag in dataset.get("tags", []):
            if ":" not in tag:
                continue
            key, value = tag.split(":")
            if key == "name":
                name = value
        print(f"Fetching dataset '{name}' with ID {dataset_id}")
        print(f"Dataset info:", pprint.pformat(dataset))
        output_path = output_dir / f"{name}.fasta"
        gi.datasets.download_dataset(
            dataset_id, file_path=str(output_path), use_default_filename=False
        )

    vp1_tree_dataset = get_datasets_from_history(
        galaxy=gi,
        history_name="TZ Red Eyes 5",
        dataset_name_prefix="IQ-TREE on data 30 and data 32: MaxLikelihood Tree",
    )
    gi.datasets.download_dataset(
        vp1_tree_dataset.get("id"),
        file_path=str(output_dir / "vp1_tree.nwk"),
        use_default_filename=False,
    )

    full_genome_tree_dataset = get_datasets_from_history(
        galaxy=gi,
        history_name="TZ Red Eyes 4",
        dataset_name_prefix="IQ-TREE on data 8 and data 93: MaxLikelihood Tree"
    )
    gi.datasets.download_dataset(
        full_genome_tree_dataset.get("id"),
        file_path=str(output_dir / "full_genome_tree.nwk"),
        use_default_filename=False
    )
