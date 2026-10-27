# DefenseFinder — genome coordinate integration

Это проверенные координаты, извлеченные из заголовков Prodigal FASTA внутри исходного архива DefenseFinder. 9 систем / 23 белковых компонента.

Файлы:
- defense_genes_with_coordinates.tsv — 23 компонента и исходные статистики моделей.
- defense_systems_with_coordinates.tsv — 9 регионов и заметки интерпретации.
- defense_system_regions.gff3 — интервалы-кандидаты, НЕ утвержденная функциональная активность.
- defensefinder_component_proteins.faa — последовательности этих 23 компонентов из запуска DefenseFinder.
- 04_defensefinder_ncbi_mapping.ipynb — проверка связи с NCBI gene/locus_tag.

Инструкция:
1. Скопируйте TSV/GFF3/FAA в ~/Desktop/bioinformatics/projects/genome_annotation/results/specialized/defensefinder/
2. Скопируйте .ipynb в ~/Desktop/bioinformatics/projects/genome_annotation/notebooks/
3. Убедитесь, что data/processed/NC_000913.3.gbk существует (сохранен в предыдущих шагах).
4. Откройте ноутбук в VS Code с kernel genome-base и выполните ячейки сверху вниз.
5. Получится файл defensefinder_ncbi_gene_mapping.tsv и список неоднозначных сопоставлений.

Важно: геномные координаты 1-based inclusive, длина референса 4,641,652 bp. Системные интервалы — min/max координат белковых компонентов, не экспериментально доказанные границы.
Нельзя трактовать MODEL wholeness и score как прямые доказательства биологической активности.
