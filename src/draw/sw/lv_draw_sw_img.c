/**
 * @file lv_draw_img.c
 *
 */

/*********************
 *      INCLUDES
 *********************/
#include "lv_draw_sw.h"
#include "../lv_img_cache.h"
#include "../../hal/lv_hal_disp.h"
#include "../../misc/lv_log.h"
#include "../../core/lv_refr.h"
#include "../../misc/lv_mem.h"
#include "../../misc/lv_math.h"

/*********************
 *      DEFINES
 *********************/
#define MAX_BUF_SIZE (uint32_t) lv_disp_get_hor_res(_lv_refr_get_disp_refreshing())

/**********************
 *      TYPEDEFS
 **********************/

/**********************
 *  STATIC PROTOTYPES
 **********************/

/**********************
 *  STATIC VARIABLES
 **********************/

/**********************
 *      MACROS
 **********************/

/**********************
 *   GLOBAL FUNCTIONS
 **********************/


/* Separate the image channels to RGB and Alpha to match LV_COLOR_DEPTH settings*/
void convert_cb(const lv_area_t * dest_area, const void * src_buf, lv_coord_t src_w, lv_coord_t src_h,
                lv_coord_t src_stride, const lv_draw_img_dsc_t * draw_dsc, lv_img_cf_t cf, lv_color_t * cbuf, lv_opa_t * abuf)
{

    const uint8_t * src_tmp8 = (const uint8_t *)src_buf;
    lv_coord_t y;
    lv_coord_t x;

    if(cf == LV_IMG_CF_RGB || cf == LV_IMG_CF_RGB_CHK) {
        uint32_t px_cnt = lv_area_get_size(dest_area);
        lv_memset_ff(abuf, px_cnt);

        src_tmp8 += (src_stride * dest_area->y1 * sizeof(lv_color_t)) + dest_area->x1 * sizeof(lv_color_t);
        uint32_t dest_w = lv_area_get_width(dest_area);
        uint32_t dest_w_byte = dest_w * sizeof(lv_color_t);

        lv_coord_t src_stride_byte = src_stride * sizeof(lv_color_t);
        lv_color_t * cbuf_tmp = cbuf;
        for(y = dest_area->y1; y <= dest_area->y2; y++) {
            lv_memcpy(cbuf_tmp, src_tmp8, dest_w_byte);
            src_tmp8 += src_stride_byte;
            cbuf_tmp += dest_w;
        }

        /*Make "holes" for with Chroma keying*/
        if(cf == LV_IMG_CF_RGB_CHK) {
            uint32_t i;
            lv_color_t chk = LV_COLOR_CHROMA_KEY;
#if LV_COLOR_DEPTH == 8
            uint8_t * cbuf_uint = (uint8_t *)cbuf;
            uint8_t chk_v = chk.full;
#elif LV_COLOR_DEPTH == 16
            uint16_t * cbuf_uint = (uint16_t *)cbuf;
            uint16_t chk_v = chk.full;
#elif LV_COLOR_DEPTH == 32
            uint32_t * cbuf_uint = (uint32_t *)cbuf;
            uint32_t chk_v = chk.full;
#endif
            for(i = 0; i < px_cnt; i++) {
                if(chk_v == cbuf_uint[i]) abuf[i] = 0x00;
            }
        }
    }
    else if(cf == LV_IMG_CF_RGBA) {
        src_tmp8 += (src_stride * dest_area->y1 * LV_IMG_PX_SIZE_ALPHA_BYTE) + dest_area->x1 * LV_IMG_PX_SIZE_ALPHA_BYTE;

        lv_coord_t src_new_line_step_px = (src_stride - lv_area_get_width(dest_area));
        lv_coord_t src_new_line_step_byte = src_new_line_step_px * LV_IMG_PX_SIZE_ALPHA_BYTE;

        lv_coord_t dest_h = lv_area_get_height(dest_area);
        lv_coord_t dest_w = lv_area_get_width(dest_area);
        for(y = 0; y < dest_h; y++) {
            for(x = 0; x < dest_w; x++) {
                abuf[x] = src_tmp8[LV_IMG_PX_SIZE_ALPHA_BYTE - 1];
#if LV_COLOR_DEPTH == 8
                cbuf[x] = *src_tmp8;
#elif LV_COLOR_DEPTH == 16
                cbuf[x].full = *src_tmp8 + ((*(src_tmp8 + 1)) << 8);
#elif LV_COLOR_DEPTH == 32
                cbuf[x] = *((lv_color_t *) src_tmp8);
                cbuf[x].ch.alpha = 0xff;
#endif
                src_tmp8 += LV_IMG_PX_SIZE_ALPHA_BYTE;

            }
            cbuf += dest_w;
            abuf += dest_w;
            src_tmp8 += src_new_line_step_byte;
        }
    }
}

static void transform_point(const lv_img_transform_dsc_t * dsc, int32_t x, int32_t y, int32_t * xs, int32_t * ys)
{
    /*Get the target point relative coordinates to the pivot*/
    int32_t xt = x - dsc->cfg.pivot_x;
    int32_t yt = y - dsc->cfg.pivot_y;

    if(dsc->cfg.zoom == LV_IMG_ZOOM_NONE) {
        /*Get the source pixel from the upscaled image*/
        *xs = ((dsc->tmp.cosma * xt - dsc->tmp.sinma * yt) >> (_LV_TRANSFORM_TRIGO_SHIFT - 8)) + dsc->tmp.pivot_x_256;
        *ys = ((dsc->tmp.sinma * xt + dsc->tmp.cosma * yt) >> (_LV_TRANSFORM_TRIGO_SHIFT - 8)) + dsc->tmp.pivot_y_256;
    }
    else if(dsc->cfg.angle == 0) {
        xt = (int32_t)((int32_t)xt * dsc->tmp.zoom_inv) >> _LV_ZOOM_INV_UPSCALE;
        yt = (int32_t)((int32_t)yt * dsc->tmp.zoom_inv) >> _LV_ZOOM_INV_UPSCALE;
        *xs = xt + dsc->tmp.pivot_x_256;
        *ys = yt + dsc->tmp.pivot_y_256;
    }
    else {
        xt = (int32_t)((int32_t)xt * dsc->tmp.zoom_inv) >> _LV_ZOOM_INV_UPSCALE;
        yt = (int32_t)((int32_t)yt * dsc->tmp.zoom_inv) >> _LV_ZOOM_INV_UPSCALE;
        *xs = ((dsc->tmp.cosma * xt - dsc->tmp.sinma * yt) >> (_LV_TRANSFORM_TRIGO_SHIFT)) + dsc->tmp.pivot_x_256;
        *ys = ((dsc->tmp.sinma * xt + dsc->tmp.cosma * yt) >> (_LV_TRANSFORM_TRIGO_SHIFT)) + dsc->tmp.pivot_y_256;
    }
}

#define CENTER_LIMIT 0x40
void tranform_cb(const lv_area_t * dest_area, const void * src_buf, lv_coord_t src_w, lv_coord_t src_h,
                 lv_coord_t src_stride, const lv_draw_img_dsc_t * draw_dsc, lv_img_cf_t cf, lv_color_t * cbuf, lv_opa_t * abuf)
{
    lv_img_transform_dsc_t trans_dsc;
    lv_memset_00(&trans_dsc, sizeof(lv_img_transform_dsc_t));
    trans_dsc.cfg.angle = draw_dsc->angle;
    trans_dsc.cfg.zoom = draw_dsc->zoom;
    trans_dsc.cfg.src = src_buf;
    trans_dsc.cfg.src_w = src_w;
    trans_dsc.cfg.src_h = src_h;
    trans_dsc.cfg.cf = cf;
    trans_dsc.cfg.pivot_x = draw_dsc->pivot.x;
    trans_dsc.cfg.pivot_y = draw_dsc->pivot.y;
    trans_dsc.cfg.color = draw_dsc->recolor;
    trans_dsc.cfg.antialias = draw_dsc->antialias;
    _lv_img_buf_transform_init(&trans_dsc);

    lv_coord_t dest_w = lv_area_get_width(dest_area);
    lv_coord_t dest_h = lv_area_get_height(dest_area);
    lv_coord_t y;
    for(y = 0; y < dest_h; y++) {
        int32_t xs1_ups, ys1_ups, xs2_ups, ys2_ups;
        transform_point(&trans_dsc, dest_area->x1, dest_area->y1 + y, &xs1_ups, &ys1_ups);
        transform_point(&trans_dsc, dest_area->x2, dest_area->y1 + y, &xs2_ups, &ys2_ups);

        int32_t xs_step = (xs2_ups - xs1_ups) / dest_w;
        int32_t ys_step = (ys2_ups - ys1_ups) / dest_w;

        lv_coord_t x;
        int32_t xs_ups = xs1_ups + xs_step / 2;     /*Init. + go the center of the pixel*/
        int32_t ys_ups = ys1_ups + ys_step / 2;
        if(trans_dsc.cfg.antialias == 0) {
            for(x = 0; x < dest_w; x++) {
                int32_t xs_int = xs_ups >> 8;
                int32_t ys_int = ys_ups >> 8;
                if(xs_int < 0 || xs_int >= src_w || ys_int < 0 || ys_int >= src_h) {
                    abuf[x] = 0;
                }
                else {
                    const uint8_t * src_tmp;
                    src_tmp = trans_dsc.cfg.src;
                    src_tmp += (ys_int * src_stride * LV_IMG_PX_SIZE_ALPHA_BYTE) + xs_int * LV_IMG_PX_SIZE_ALPHA_BYTE;
                    abuf[x] = src_tmp[LV_IMG_PX_SIZE_ALPHA_BYTE - 1];
                    cbuf[x].full = src_tmp[0] + (src_tmp[1] << 8);
                }
                xs_ups += xs_step;
                ys_ups += ys_step;
            }
        }
        else {
            for(x = 0; x < dest_w; x++) {
                int32_t xs_int = xs_ups >> 8;
                int32_t ys_int = ys_ups >> 8;
                int32_t xs_fract = xs_ups & 0xFF;
                int32_t ys_fract = ys_ups & 0xFF;

                const uint8_t * src_tmp;
                src_tmp = trans_dsc.cfg.src;
                src_tmp += (ys_int * src_stride * LV_IMG_PX_SIZE_ALPHA_BYTE) + xs_int * LV_IMG_PX_SIZE_ALPHA_BYTE;

                int32_t x_next;
                int32_t y_next;
                const uint8_t * px_base = NULL;
                if(xs_int >= 0 && xs_int < src_w && ys_int >= 0 && ys_int < src_h) {
                    /*We hit the ~center of a pixel, get it directly*/
                    if(xs_fract > 0x80 - CENTER_LIMIT &&
                       xs_fract < 0x80 + CENTER_LIMIT &&
                       ys_fract > 0x80 - CENTER_LIMIT &&
                       ys_fract < 0x80 + CENTER_LIMIT) {
                        abuf[x] = src_tmp[LV_IMG_PX_SIZE_ALPHA_BYTE - 1];
                        cbuf[x].full = src_tmp[0] + (src_tmp[1] << 8);
                        ys_ups += ys_step;
                        xs_ups += xs_step;
                        continue;
                    }
                    /*Edge of the pixel, make some expensive filtering*/
                    else {
                        px_base = src_tmp;
                    }
                }
                /*Fully out of the image*/
                else if(xs_int < -1 || xs_int > src_w || ys_int < -1 || ys_int > src_h) {
                    abuf[x] = 0x00;
                    ys_ups += ys_step;
                    xs_ups += xs_step;
                    continue;
                }

                //            abuf[x] = 0xff;
                //            cbuf[x].full = 0xf000;
                //            ys_ups += ys_step;
                //            xs_ups += xs_step;
                //            continue;


                /*Get the direction the hor and ver neighbor
                 *`fract` will be in range of 0x00..0xFF and `next` (+/-1) indicates the direction*/
                if(xs_fract < 0x80) {
                    x_next = -1;
                    xs_fract = (0x7F - xs_fract) * 2;
                }
                else {
                    x_next = 1;
                    xs_fract = (xs_fract - 0x80) * 2;
                }
                if(ys_fract < 0x80) {
                    y_next = -1;
                    ys_fract = (0x7F - ys_fract) * 2;
                }
                else {
                    y_next = 1;
                    ys_fract = (ys_fract - 0x80) * 2;
                }

                /*Get the horizontal neighbor*/
                const uint8_t * px_hor = NULL;
                if(ys_int >= 0 && ys_int < src_h && xs_int + x_next >= 0 && xs_int + x_next < src_w) {
                    px_hor = src_tmp + x_next * LV_IMG_PX_SIZE_ALPHA_BYTE;
                }

                /*Get the vertical neighbor*/
                const uint8_t * px_ver = NULL;
                if(ys_int + y_next >= 0 && ys_int + y_next < src_h && xs_int >= 0 && xs_int < src_w) {
                    px_ver = src_tmp + y_next * src_stride * LV_IMG_PX_SIZE_ALPHA_BYTE;
                }

                if(px_base == NULL || px_hor == NULL || px_ver == NULL) {
                    //                    ys_ups += ys_step;
                    //                    xs_ups += xs_step;
                    //                    continue;
                    /*If all are out of the image just set alpha to 0x00*/
                    //                    printf("---\n");
                    if(px_base == NULL && px_hor == NULL && px_ver == NULL) {
                        abuf[x] = 0x00;
                        ys_ups += ys_step;
                        xs_ups += xs_step;
                        continue;
                    }
                    /*More complex cases, partially on the image*/
                    else if(px_base == NULL && px_hor == NULL) {
                        /*On the top/bottom edge*/
                        px_base = px_ver;
                        px_hor = px_ver;
                    }
                    else if(px_base == NULL && px_ver == NULL) {
                        /*On the left/right edge*/
                        px_base = px_hor;
                        px_ver = px_hor;
                    }
                    else if(px_hor == NULL && px_ver == NULL) {
                        /*On the corners*/
                        px_hor = px_base;
                        px_ver = px_base;
                    }

                    //                printf("%d, %d, %d, %d\n", xs_int, ys_int, x_next, y_next);

                    /*Now only hor or ver can be NULL
                     *It's not possible that base is NULL is hor and ver was on the image*/
                    if(px_ver == NULL) {
                        px_ver = px_base;

                    }
                    else if(px_hor == NULL) {
                        px_hor = px_base;
                    }
                }
                else {
                    //                    printf("-\n");
                }
                //                    printf("\n");
                //
                lv_opa_t a_base = px_base[LV_IMG_PX_SIZE_ALPHA_BYTE - 1];
                lv_opa_t a_ver = px_ver[LV_IMG_PX_SIZE_ALPHA_BYTE - 1];
                lv_opa_t a_hor = px_hor[LV_IMG_PX_SIZE_ALPHA_BYTE - 1];

                if(a_ver != a_base) a_ver = ((a_ver * ys_fract) + (a_base * (0xFF - ys_fract))) >> 8;
                if(a_hor != a_base) a_hor = ((a_hor * xs_fract) + (a_base * (0xFF - xs_fract))) >> 8;
                abuf[x] = (a_ver + a_hor) >> 1;

                if(abuf[x]) {
                    lv_color_t c_base;
                    lv_color_t c_ver;
                    lv_color_t c_hor;
                    c_base.full = px_base[0] + (px_base[1] << 8);
                    c_ver.full = px_ver[0] + (px_ver[1] << 8);
                    c_hor.full = px_hor[0] + (px_hor[1] << 8);
                    c_ver = lv_color_mix(c_ver, c_base, ys_fract);
                    c_hor = lv_color_mix(c_hor, c_base, xs_fract);
                    cbuf[x] = c_base; //lv_color_mix(c_hor, c_ver, LV_OPA_50);
                }

                xs_ups += xs_step;
                ys_ups += ys_step;
            }
        }
        cbuf += dest_w;
        abuf += dest_w;
    }
}


LV_ATTRIBUTE_FAST_MEM void lv_draw_sw_img_decoded(struct _lv_draw_ctx_t * draw_ctx, const lv_draw_img_dsc_t * draw_dsc,
                                                  const lv_area_t * coords, const uint8_t * src_buf, lv_img_cf_t cf)
{
    /*Use the clip area as draw area*/
    lv_area_t draw_area;
    lv_area_copy(&draw_area, draw_ctx->clip_area);

    bool mask_any = lv_draw_mask_is_any(&draw_area);
    bool transform = draw_dsc->angle != 0 || draw_dsc->zoom != LV_IMG_ZOOM_NONE ? true : false;


    lv_area_t blend_area;
    lv_draw_sw_blend_dsc_t blend_dsc;

    lv_memset_00(&blend_dsc, sizeof(lv_draw_sw_blend_dsc_t));
    blend_dsc.opa = draw_dsc->opa;
    blend_dsc.blend_mode = draw_dsc->blend_mode;
    blend_dsc.blend_area = &blend_area;

    /*The simplest case just copy the pixels into the draw_buf*/
    if(!mask_any && !transform && cf == LV_IMG_CF_RGB && draw_dsc->recolor_opa == LV_OPA_TRANSP) {
        blend_dsc.src_buf = (const lv_color_t *)src_buf;

        blend_dsc.blend_area = coords;
        lv_draw_sw_blend(draw_ctx, &blend_dsc);
    }
    /*In the other cases every pixel need to be checked one-by-one*/
    else {
        blend_area.x1 = draw_ctx->clip_area->x1;
        blend_area.x2 = draw_ctx->clip_area->x2;
        blend_area.y1 = draw_ctx->clip_area->y1;
        blend_area.y2 = draw_ctx->clip_area->y2;

        lv_coord_t src_w = lv_area_get_width(coords);
        lv_coord_t src_h = lv_area_get_height(coords);
        lv_coord_t blend_h = lv_area_get_height(&blend_area);
        lv_coord_t blend_w = lv_area_get_width(&blend_area);

        uint32_t max_buf_size = MAX_BUF_SIZE;
        uint32_t blend_size = lv_area_get_size(&blend_area);
        uint32_t buf_h;
        uint32_t buf_w = blend_w;
        if(blend_size <= max_buf_size) {
            buf_h = blend_h;
        }
        else {
            /*Round to full lines*/
            buf_h = max_buf_size / blend_w;
        }

        /*Create buffers and masks*/
        uint32_t buf_size = buf_w * buf_h;

        lv_color_t * rgb_buf = lv_mem_buf_get(buf_size * sizeof(lv_color_t));
        lv_opa_t * mask_buf = lv_mem_buf_get(buf_size);
        blend_dsc.mask_buf = mask_buf;
        blend_dsc.mask_area = &blend_area;
        blend_dsc.mask_res = LV_DRAW_MASK_RES_CHANGED;
        blend_dsc.src_buf = rgb_buf;
        lv_coord_t y_last = blend_area.y2;
        blend_area.y2 = blend_area.y1 + buf_h - 1;

        bool transform = draw_dsc->angle != 0 || draw_dsc->zoom != LV_IMG_ZOOM_NONE ? true : false;

        lv_draw_mask_res_t mask_res_def = (cf != LV_IMG_CF_RGB || draw_dsc->angle || draw_dsc->zoom != LV_IMG_ZOOM_NONE) ?
                                          LV_DRAW_MASK_RES_CHANGED : LV_DRAW_MASK_RES_FULL_COVER;
        blend_dsc.mask_res = mask_res_def;

        while(blend_area.y1 <= y_last) {
            /*Apply transformations if any or separate the channels*/
            lv_area_t transform_area;
            lv_area_copy(&transform_area, &blend_area);
            lv_area_move(&transform_area, -coords->x1, -coords->y1);
            if(transform) {
                tranform_cb(&transform_area, src_buf, src_w, src_h, src_w, draw_dsc, cf, rgb_buf, mask_buf);
            }
            else {
                convert_cb(&transform_area, src_buf, src_w, src_h, src_w, draw_dsc, cf, rgb_buf, mask_buf);
            }

            /*Apply recolor*/
            if(draw_dsc->recolor_opa > LV_OPA_MIN) {
                uint16_t premult_v[3];
                lv_opa_t recolor_opa = draw_dsc->recolor_opa;
                lv_color_t recolor = draw_dsc->recolor;
                lv_color_premult(recolor, recolor_opa, premult_v);
                uint32_t i;
                for(i = 0; i < buf_size; i++) {
                    rgb_buf[i] = lv_color_mix_premult(premult_v, rgb_buf[i], recolor_opa);
                }
            }

            /*Apply the masks if any*/
            if(mask_any) {
                lv_coord_t y;
                lv_opa_t * mask_buf_tmp = mask_buf;
                for(y = blend_area.y1; y <= blend_area.y2; y++) {
                    lv_draw_mask_res_t mask_res_line;
                    mask_res_line = lv_draw_mask_apply(mask_buf_tmp, blend_area.x1, y, blend_w);

                    if(mask_res_line == LV_DRAW_MASK_RES_TRANSP) {
                        lv_memset_00(mask_buf_tmp, blend_w);
                        blend_dsc.mask_res = LV_DRAW_MASK_RES_CHANGED;
                    }
                    else if(mask_res_line == LV_DRAW_MASK_RES_CHANGED) {
                        blend_dsc.mask_res = LV_DRAW_MASK_RES_CHANGED;
                    }
                    mask_buf_tmp += blend_w;
                }
            }

            /*Blend*/
            lv_draw_sw_blend(draw_ctx, &blend_dsc);

            /*Go the the next lines*/
            blend_area.y1 = blend_area.y2 + 1;
            blend_area.y2 = blend_area.y1 + buf_h - 1;
            if(blend_area.y2 > y_last) blend_area.y2 = y_last;
        }

        lv_mem_buf_release(mask_buf);
        lv_mem_buf_release(rgb_buf);
    }
}

/**********************
 *   STATIC FUNCTIONS
 **********************/
