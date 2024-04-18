import shutil


def copy_changelog_and_license(*args, **kwargs):
    shutil.copy("CHANGELOG.md", "docs/changelog.md")
    shutil.copy("LICENSE.md", "docs/license.md")
    shutil.copy("config/config.yaml", "docs/includes/config.yaml")
    shutil.copy("config/multiqc_config_DNA.yaml", "docs/includes/multiqc_config.yaml")
    shutil.copy("config/resources.yaml", "docs/includes/resources.yaml")
    shutil.copy("images/dag.svg", "docs/includes/images/dag.svg")