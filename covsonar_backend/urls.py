from django.urls import path, include
from rest_framework import routers, serializers, viewsets
from rest_api import viewsets
from django.contrib import admin
from . import settings

router = routers.DefaultRouter()
router.register(r"variants", viewsets.VariantViewSet)
router.register(r"genes", viewsets.GenesViewSet)
router.register(r"snp1", viewsets.SNP1ViewSet)
router.register(r"mutation_signature", viewsets.MutationSignatureViewSet)
router.register(r"samples", viewsets.SampleViewSet)
router.register(r"sample_genomes", viewsets.SampleGenomeViewSet)
router.register(r"references", viewsets.ReferenceViewSet)
router.register(r"elements", viewsets.ElementViewSet)
router.register(r"property", viewsets.PropertyViewSet)
router.register(r"aa_mutations", viewsets.AAMutationViewSet)

urlpatterns = [
    path("admin/", admin.site.urls),
    path("", include(router.urls)),
]
if settings.DEBUG:
    import debug_toolbar

    urlpatterns += [
        path('__debug__/', include(debug_toolbar.urls)),
    ]