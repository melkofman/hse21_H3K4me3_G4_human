# hse21_H3K4me3_G4_human
## HSE bioinf project
### Отчет по работе над проектом. 
### Студент: Кофман М.О., группа 3
#### Начало работы 
С помощью команды wget были получены необходимые файлы:

`wget https://www.encodeproject.org/files/ENCFF023LTU/@@download/ENCFF023LTU.bed.gz`

`wget https://www.encodeproject.org/files/ENCFF432EMI/@@download/ENCFF432EMI.bed.gz`

Оставили первые пять столбцов:

`zcat ENCFF023LTU.bed.gz | cut -f1-5 > H3K4me3_GM12878.ENCFFo23LTU.hg38.bed`

`zcat ENCFF432EMI.bed.gz | cut -f1-5 > H3K4me3_GM12878.ENCFF432EMI.hg38.bed`

С помощью liftOver преобразовали координаты из 38 в 19 версию генома:

`liftOver H3K4me3_GM12878.ENCFFo23LTU.hg38.bed hg38ToHg19.over.chain.gz H3K4me3_GM12878.ENCFF023LTU.hg19.bed H3K4me3_GM12878.ENCFF023LTU.unmapped.bed`

`liftOver H3K4me3_GM12878.ENCFF432EMI.hg38.bed hg38ToHg19.over.chain.gz H3K4me3_GM12878.ENCFF432EMI.hg19.bed H3K4me3_GM12878.ENCFF432EMI.unmapped.bed`

Затем скачали на локальный ПК данные с сервера. 

#### Распределение длин участков.

[Скрипт](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/src/len_hist.R) позволяет получить диаграммы распределения длин участков. 

![alt text](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/images/len_hist.H3K4me3_GM12878.ENCFF023LTU.hg19.png "H3K4me3_GM12878.ENCFF023.hg19 до фильтрации")

![alt text](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/images/len_hist.H3K4me3_GM12878.ENCFFo23LTU.hg38.png "H3K4me3_GM12878.ENCFFo23LTU.hg38 до фильтрации")

![alt text](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/images/len_hist.H3K4me3_GM12878.ENCFF432EMI.hg19.png "H3K4me3_GM12878.ENCFF432EMI.hg19 до фильтрации")

![alt text](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/images/len_hist.H3K4me3_GM12878.ENCFF432EMI.hg38.png "H3K4me3_GM12878.ENCFF432EMI.hg38 до фильтрации")

Теперь отфильтруем, уберем сильно длинные участки (>50000) с помощью [скрипта](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/src/filter_peaks.R).
Опять построим распределение длин участков и посмотрим, что изменилось: 

![alt text](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/images/len_hist.H3K4me3_GM12878.ENCFF023LTU.hg19.filtered.png "H3K4me3_GM12878.ENCFF023LTU.hg19 после фильтрации")

![alt text](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/images/len_hist.H3K4me3_GM12878.ENCFFo23LTU.hg38.filtered.png "H3K4me3_GM12878.ENCFF023LTU.hg38 после фильтрации")

![alt text](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/images/len_hist.H3K4me3_GM12878.ENCFF432EMI.hg19.filtered.png "H3K4me3_GM12878.ENCFF432EMI.hg19 после фильтрации")

![alt text](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/images/len_hist.H3K4me3_GM12878.ENCFF432EMI.hg38.filtered.png "H3K4me3_GM12878.ENCFF432EMI.hg38 после фильтрации")

Как видим, сильных изменений не произошло, было отрезано всего лишь одно чтение. 
Для дальнейшей работы будем использовать отфильтрованные файлы 19 версии. 

Визуализируем в геномном браузере: 

`track visibility=dense name="ENCFF023LTU"  description="H3K4me3_GM12878.ENCFF023LTU.hg19.filtered.bed"
https://raw.githubusercontent.com/melkofman/bioinf_project_G4_H3K4me3/mokofman/data/H3K4me3_GM12878.ENCFF023LTU.hg19.bed`

`track visibility=dense name="ENCFF432EMI" description="H3K4me3_GM12878.ENCFF432EMI.hg19.filtered.bed"
https://raw.githubusercontent.com/melkofman/bioinf_project_G4_H3K4me3/mokofman/data/filtered/H3K4me3_GM12878.ENCFF432EMI.hg19.filtered.bed`

Как видно из визуализации, есть участки аннотации генов, которые попадают в промоторные области. 

#### Куда попадают участки по аннотации генов. 
Будем использовать библиотеку ChIPseeker и [скрипт](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/src/peakAnno.R)

Результаты: 

![alt text](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/images/chip_seeker.H3K4me3_GM12878.ENCFF023LTU.hg19.filtered.plotAnnoPie.png "H3K4me3_GM12878.ENCFF023LTU.hg19")

Как видно по диаграмме для H3K4me3_GM12878.ENCFF023LTU.hg19, бОльшая часть попадает в промоторы. 

![alt text](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/images/chip_seeker.H3K4me3_GM12878.ENCFF432EMI.hg19.filtered.plotAnnoPie.png "H3K4me3_GM12878.ENCFF432EMI.hg19")

Для H3K4me3_GM12878.ENCFF432EMI.hg19 картина другая, бОльшая часть попадает в интроны и экзоны. 

Ниже приведены также распределения по хромосомам. Но так-как на каждую хромосому попадает достаточно большое количество чтений, то картинки не являются сильно информативными. 

![alt text](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/images/chip_seeker.H3K4me3_GM12878.ENCFF023LTU.hg19.filtered.covplot.pdf "H3K4me3_GM12878.ENCFF023LTU.hg19 по хромосомам")
![alt text](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/images/chip_seeker.H3K4me3_GM12878.ENCFF432EMI.hg19.filtered.covplot.pdf "H3K4me3_GM12878.ENCFF432EMI.hg19 по хромосомам")


#### Работа со вторичной структурой ДНК.
Вторичная структура G4 находится в [файле](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/data/G4.bed)

Визуализируем в геномном браузере: 

`track visibility=dense name="G4"  color=0,200,0  description="G4"
https://raw.githubusercontent.com/melkofman/bioinf_project_G4_H3K4me3/mokofman/data/G4.bed
`
![alt text](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/images/G4_genBr.png "визуализация G4 в геномном браузере")

Построим распределение длин G4:
![alt text](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/images/len_hist.G4.png " распределение длины")


Теперь посмотрим, куда попадают участки по аннотации генов, используя тот же скрипт, что и выше:

![alt text](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/images/chip_seeker.G4.plotAnnoPie.png "участки по аннотации")

#### Анализ пересечений гистоновой метки и структуры ДНК
Объединим пики

`cat  *hg19.filtered.bed  |   sort -k1,1 -k2,2n   |   bedtools merge   >  H3K4me3_GM12878.merge.hg19.bed`


 Визуализируем это объединение в геномном браузере: 
 
 `track visibility=dense name="ChIP_merge"  color=50,50,200   description=" H3K4me3_GM12878.merge.hg19.bed"
https://raw.githubusercontent.com/melkofman/bioinf_project_G4_H3K4me3/mokofman/data/filtered/H3K4me3_GM12878.merge.hg19.bed
`

 ![alt text](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/images/chipmerge.png "объединенные пики")
 
 Пересекаем пики гистоновой метки со вторичной стуктурой ДНК: 
 
 
 `bedtools intersect  -a /Users/melanie/Desktop/bioinf/project/gitpr/bioinf_project_G4_H3K4me3/data/G4.bed   -b  /Users/melanie/Desktop/bioinf/project/gitpr/bioinf_project_G4_H3K4me3/data/filtered/H3K4me3_GM12878.merge.hg19.bed >  H3K4me3_GM12878_intersect_G4.bed
 `
 
 
 
 И визуализируем это: 
 
 `track visibility=dense name="intersect_with_G4"  color=200,0,0  description=" H3K4me3_GM12878_intersect_G4.bed"
https://raw.githubusercontent.com/melkofman/bioinf_project_G4_H3K4me3/mokofman/data/H3K4me3_GM12878_intersect_G4.bed
`


 ![alt text](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/images/int_g_b.png "пересечение")
 Координаты: chr1:762,047-762,071
 
 ![alt text](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/images/intersect_genbr.png "пересечение")
 
 
 Красным отмечено пересечение, которое мы получили. Синим - объединенные пики. Зеленым - вторичная структура ДНК. 
 Видно, что красная полоса соответствует пересечению синей и зеленой. 
 
 [Сессия ген браузера](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/data/gen_browser_session) 
#### Ассоциация пиков с генами. 
  [Скрипт](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/src/ChIpPeakAnno.R) позволяет получить ассоциацию пиков с генами.  
 
 Результаты: [гены](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/data/H3K4me3_GM12878_intersect_G4.genes_uniq.txt), [уникальный набор генов](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/data/H3K4me3_GM12878_intersect_G4.genes.txt)
 
#### GO анализ
 С помощью [Panther](http://pantherdb.org/) проведем GO анализ. 
 
 Зададим следующие настройки для анализа: 
 
 ![alt text](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/images/pan_setup.png "настройки для go анализа")
 
 Результаты: 
 
 ![alt text](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/images/resultGO.png "результаты go анализа")
 
 Всего генов 20595. В исходном списке было 943 гена, из которых програма распознала 913. 
 
 [Результаты в виде текстового документа](https://github.com/melkofman/bioinf_project_G4_H3K4me3/blob/mokofman/data/analysis.txt)
 
 
 
 
 
 
