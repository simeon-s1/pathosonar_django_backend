
# HC:
# Frameshifts []<VariantLabel>:
#   v.frameshift = 1
# Genomic_profiles []<VariantLabel>:
#    element_type: nt
#    variant_alt != 'N'
# proteomic_profiles []<VariantLabel>:
# -  element_type: cds
#    variant_alt != 'X'
# 
# INPUT:
# Properties
# Variant - start, ende etc..





SELECT fs.sample_id AS 'sample.id',
       fs.name AS 'sample.name',
       fs.seqhash,       
       GROUP_CONCAT(CASE
                        WHEN v.element_id = 1
                             AND v.frameshift = 1 THEN v.label
                    END, ' ') AS frameshifts,
       GROUP_CONCAT(CASE
                        WHEN v.element_id = 1
                             AND v.alt != 'N'
                             AND v.alt != '.' THEN v.label
                    END, ' ') AS genomic_profiles,
       GROUP_CONCAT(CASE
                        WHEN v.element_id IN (12, 13, 14, 15, 16, 17, 18, 19, 20, 21)
                             AND v.alt != 'X' THEN e.symbol || ':' || v.label
                    END, ' ') AS proteomic_profiles
FROM
  (SELECT sub.sample_id,
          sub.name,
          sub.seqhash
   FROM
     (SELECT s.id AS 'sample_id',
             s.name,
             s.seqhash,
             CASE
                 WHEN sample2property.property_id = 3
                      AND sample2property.value_date BETWEEN '2020-01-01' AND '2023-12-31' THEN 1
                 ELSE 0
             END AS property_1,
             CASE
                 WHEN molecule.standard = 1
                      AND element.type = 'cds'
                      AND element.symbol = 'S'
                      AND variant.start = 500
                      AND variant.end = 501
                      AND variant.ref = 'N'
                      AND variant.alt = 'Y' THEN 1
                 ELSE 0
             END AS mutation_1
      FROM sample s
      JOIN sample2property ON s.id = sample2property.sample_id
      JOIN alignment ON s.seqhash = alignment.seqhash
      JOIN alignment2variant ON alignment.id = alignment2variant.alignment_id
      JOIN variant ON alignment2variant.variant_id = variant.id
      JOIN element ON variant.element_id = element.id
      JOIN molecule ON element.molecule_id = molecule.id
      WHERE property_1 = 1
        AND mutation_1 = 1) AS sub) AS fs
LEFT JOIN alignment a ON fs.seqhash = a.seqhash
LEFT JOIN alignment2variant a2v ON a.id = a2v.alignment_id
LEFT JOIN variant v ON a2v.variant_id = v.id
LEFT JOIN element e ON v.element_id = e.id
GROUP BY fs.sample_id
"""

