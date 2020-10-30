#include <openssl/bn.h>
#include "hesse.h"
#include <stdio.h>
#include <string.h>

int main()
{
    struct par param = {NULL, NULL, NULL, NULL, NULL, NULL};
    par_init(&param);

    //вывожу наборы параметров
    printf("Набор параметров:\n");
    printf("p=%s\n", BN_bn2dec(param.p));
    printf("a=%s\n", BN_bn2dec(param.a));
    printf("u=%s\n", BN_bn2dec(param.u));
    printf("v=%s\n", BN_bn2dec(param.v));

    //параметры Скрученного Гессе

    printf("Параметры Скрученного Гессе:\n");
    struct twisted_hesse curve={NULL, NULL, NULL, NULL, NULL, NULL};
    twisted_hesse_init(&curve, &param);
    printf("a =%s\n", BN_bn2dec(curve.a));
    printf("d =%s\n", BN_bn2dec(curve.d));
    printf("X_base=%s\n", BN_bn2dec(curve.X));
    printf("Y_base=%s\n", BN_bn2dec(curve.Y));
    printf("Z_base=%s\n\n", BN_bn2dec(curve.Z));

    printf("Нейтральный элемент:\n");
    struct point O = {NULL, NULL, NULL};
    point_init(&O ,"0", "-1", "1");
    printf("в проективных координатах:\n");
    print_in_projective(&O);
    printf("в афинных:\n");
    print_in_affine(&O);

    if(aff_point_check(&O, &curve)){
        printf("Точка O находится на кривой\n\n");
    }
    else{
        printf("Точка O не находится на кривой\n\n");
    }

    printf("Тест 1:\n");
    printf("Лежит ли точка P=(2:3:4) на кривой\n");
    printf("P в аффинных координатах:\n");
    struct point P2 = {NULL, NULL, NULL};
    point_init(&P2, "2", "3", "4");
    print_in_affine(&P2);
    if(aff_point_check(&P2,&curve)){
        printf("Точка P находится на кривой\n\n");
    }
    else{
        printf("Точка P не находится на кривой\n\n");
    }

    printf("Тест 2: Проверка [q]P = O\n");
    struct point P_s = {NULL, NULL, NULL};
    point_init(&P_s,"0","-1","1");
    struct point res_point = {NULL, NULL, NULL};
    point_init(&res_point, "0", "-1" ,"1");
    printf("Нейтральный элемент в аффинных координатах:\n");
    print_in_affine(&O);
    printf("qP в аффинных координатах:\n");
    print_in_affine(&res_point);


    printf("Тест 3:\n");
    printf("проверим, что [q+1]P = P и [q-1] = -P\n");
    BIGNUM* degree = BN_new();
    BN_dec2bn(&degree, "1");// degree=1;
    BN_add (degree, degree, param.q);//degree=q+1
    crat_find(&res_point, &P_s, &curve, degree);
    point_init(&res_point, "0", "-1" ,"1");
    printf("[q+1]P:\n");
    print_in_affine(&res_point);
    printf("P:\n");
    print_in_affine(&P_s);
    if(!(is_point_equal(&res_point, &P_s,&curve))){
        printf("[q+1]P равно P\n\n");
    }
    else{
        printf("[q+1]P не равно P\n\n");
    }

    BN_dec2bn (&degree, "1");                      // degree = 1
    BN_sub (degree, param.q, degree);                // degree = q-1
    crat_find(&res_point, &P_s, &curve,degree);
    printf("[q-1]P:\n");
    print_in_affine(&res_point);
    printf("-P:\n");
    struct point reverse = {NULL, NULL, NULL};
    point_init(&reverse, "3", "6", "8");
    reverse_point(&reverse, &P_s, &param);
    print_in_affine(&reverse);
    if(!(is_point_equal(&reverse, &reverse, &curve))){
        printf("[q-1]P не равно -P\n\n");
    }
    else{
        printf("[q+1]P равно P\n\n");
    }

    printf("Тест 4:\n");
    printf("[k1]P + [k2]P = [k1 + k2]P\n");
    BIGNUM* k1 = BN_new();
    BIGNUM* k2 = BN_new();
    printf("Генерирую k1 и k2\n");
    BIGNUM* maxrand = BN_new();
    BN_dec2bn(&maxrand, "100000000000000000");
    BN_rand_range(k1, maxrand);
    BN_rand_range(k2, maxrand);
    printf("k1 = %s\n", BN_bn2dec(k1));
    printf("k2 = %s\n", BN_bn2dec(k2));
    BN_add(param.q,k1,k2);                          //k=k1+k2
    printf("k = k1 + k2 = %s\n", BN_bn2dec(param.q));
    struct point res1={NULL, NULL, NULL};
    struct point res2={NULL, NULL, NULL};
    struct point res3={NULL, NULL, NULL};
    point_init(&res1,"0", "1", "1");
    point_init(&res2,"0", "1", "1");
    point_init(&res3,"0", "1", "1");

    add_points(&res1, &res2, &res_point, &curve);

    if(!(is_point_equal(&res_point, &res3, &curve))){
        printf("[k1]P + [k2]P равно [k1 + k2]P\n\n");
    }
    else{
        printf("[k1]P + [k2]P не равно [k1 + k2]P\n\n");
    }
    if(aff_point_check(&res3, &curve)){
        printf("точка [k]P находится на кривой\n\n");
    }
    else{
        printf("точка [k]P не находится на кривой\n\n");
    }


}